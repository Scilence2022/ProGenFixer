#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h> 
#include <math.h>

#include "ketopt.h" // command-line argument parser
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop


#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
#define KC_BITS 10
#define KC_MAX ((1<<KC_BITS) - 1)
#define kc_c4_eq(a, b) ((a)>>KC_BITS == (b)>>KC_BITS) // lower 10 bits for counts; higher bits for k-mer
#define kc_c4_hash(a) ((a)>>KC_BITS)

#define Max_Path_Num 1000
//#define Max_Path_Len 1000
#define Max_Cov_Ratio 3

#define DEBUG 0

// Define a macro for debug output that is only enabled when the DEBUG flag is set
#if DEBUG
#define debug_print(fmt, ...) \
    fprintf(stdout, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#else
#define debug_print(fmt, ...) \
    do {} while (0)
#endif

KHASHL_SET_INIT(, kc_c4_t, kc_c4, uint64_t, kc_c4_hash, kc_c4_eq)

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

const unsigned char nt4_seq_table[4] = { // translate 0123 to ACGT
        'A', 'C', 'G', 'T'};


static inline uint64_t hash64(uint64_t key, uint64_t mask) // invertible integer hash function
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

// The inversion of hash64(). Modified from <https://naml.us/blog/tag/invertible>
static inline uint64_t hash64i(uint64_t key, uint64_t mask) //
{
    uint64_t tmp;
    // Invert key = key + (key << 31)
    tmp = (key - (key << 31)); 	key = (key - (tmp << 31)) & mask;
    // Invert key = key ^ (key >> 28)
    tmp = key ^ key >> 28; 	key = key ^ tmp >> 28;
    // Invert key *= 21
    key = (key * 14933078535860113213ull) & mask;
    // Invert key = key ^ (key >> 14)
    tmp = key ^ key >> 14; 	tmp = key ^ tmp >> 14; 	tmp = key ^ tmp >> 14; 	key = key ^ tmp >> 14;
    // Invert key *= 265
    key = (key * 15244667743933553977ull) & mask;
    // Invert key = key ^ (key >> 24)
    tmp = key ^ key >> 24; 	key = key ^ tmp >> 24;
    // Invert key = (~key) + (key << 21)
    tmp = ~key; 	tmp = ~(key - (tmp << 21)); tmp = ~(key - (tmp << 21)); key = ~(key - (tmp << 21)) & mask;
    return key;
}


unsigned char* uint64_acgt(uint64_t key, unsigned char* seq, unsigned char km_len){ //decode uint_64_t kmer  to actg
    int p = km_len-1;
    //MALLOC(seq, km_len+1);
    //seq[km_len+1] = '\0';
    while (p >= 0){
        int n = key % 4;
        key = key >> 2;
        *(seq + p) = nt4_seq_table[n];
        p = p - 1;
    }
    return seq;
}

unsigned char* uint64_int8(uint64_t key, unsigned char* seq){
    //decode uint_64_t kmer  to  char

    uint64_t mask = (1ULL<<8)-1;
    int i =0;
    for(i =0; i< 8; i++){
        *(seq+7-i) = key & mask;
        mask = mask <<8;
    }
    return seq;
}



uint64_t comp_rev2(uint64_t x, unsigned char km_len){
    uint64_t a, b=0ULL, mask = (1ULL<<km_len*2) - 1;

    a = ~x; a = a & mask;
    while((a & mask) > 0 ){
        b = b<<2;
        b = b + (a & 3ULL);
        a = a>>2;
    }
//    b = ~b;
    b = b & mask;
    return b;
}


uint64_t comp_rev(uint64_t x, unsigned char km_len){
//    if (c < 4) { // not an "N" base
//        x[0] = (x[0] << 2 | c) & mask;                  // forward strand
//        x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
//    uint64_t y=0;// mask = (1ULL << km_len * 2) - 1;
    uint64_t y=0;
    int i, c;
    for(i=0;i<km_len;i++){
        y = y << 2;
        c = x & 3ULL;
        y = y | (3-c);
        x = x >> 2;
    }
    return y;
}


uint64_t min_hash_key(uint64_t x, unsigned char km_len){
    uint64_t y = comp_rev(x, km_len);
    return (x < y) ? x:y;
}


uint64_t actgkmer_uint64(unsigned char* seq, unsigned char km_len){
//    Just for testing
//// Function verified 20210708 Lifu Song
    int i;
    uint64_t x, mask = (1ULL<<km_len*2) - 1; //shift = (km_len - 1) * 2
    for (i = 0, x = 0; i < km_len; ++i) {
        //fprintf(stdout, "%c ", seq[i]);
        int c = seq_nt4_table[(uint8_t)seq[i]];
        //fprintf(stdout, "%d ", c);
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & mask;                  // forward strand
//            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
        } else x = 0; // if there is an "N", restart
    }

    //fprintf(stdout, "y %ld\n", x);
    return x;
}


uint64_t actgkmer_hashkey(unsigned char* seq, unsigned char km_len){
//// Function verified 20210708 Lifu Song
    int i;
    uint64_t x[2], y,  shift = (km_len - 1) * 2, mask = (1ULL<<km_len*2) - 1;
    for (i = 0, x[0] = x[1] = 0; i < km_len; ++i) {
        //fprintf(stdout, "%c ", seq[i]);
        int c = seq_nt4_table[(uint8_t)seq[i]];
        //fprintf(stdout, "%d ", c);
        if (c < 4) { // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
        } else x[0] = x[1] = 0; // if there is an "N", restart
    }
    y = x[0] < x[1]? x[0] : x[1];
    //fprintf(stdout, "y %ld\n", y);
    return y;
}

unsigned long uint8_actg(unsigned char bts){
    unsigned int i=0;
    unsigned long four_base_seq=0;

    for(i =0 ; i < 4; ++i){
        int m = bts & 3;
        bts = bts >>2;
        four_base_seq =four_base_seq >>8;
        four_base_seq = four_base_seq + ((unsigned long)nt4_seq_table[m]<< 24);
    }
    return four_base_seq;
}


unsigned char actg_uint8(unsigned long actg_seq){
    unsigned int i=0;
    unsigned char uint8_value=0;

    for(i =0 ; i < 4; ++i){
        unsigned int m = actg_seq & 255;
        actg_seq = actg_seq >> 8;
        uint8_value = uint8_value >> 2;
        uint8_value = uint8_value + ((unsigned char)seq_nt4_table[m]<<6);
    }
    return uint8_value;
}

void actg_intseq(unsigned char* actg_seq, unsigned char* int_seq, int8_t actg_seq_len){

    int p = 0;
    unsigned long n, n1, n2, n3, n4;
    while(p < actg_seq_len-3){
        n1 = actg_seq[p];
        n2 = actg_seq[p+1];
        n3 = actg_seq[p+2];
        n4 = actg_seq[p+3];
        n = (n1<<24) + (n2<<16) + (n3<<8) + n4;
        int_seq[p>>2] = actg_uint8(n) ;
        p = p + 4;
    }
    int_seq[p + 1] = '\0';
}


void intseq_actg(unsigned char* int_seq, unsigned char* actg_seq, int8_t int_seq_len){
//// int seq to DNA seq
//// Function not tested yet
//// one char(byte) to four bases
    int p = 0;
    unsigned long n, mask;

    unsigned char n1, n2, n3, n4;
    while(p < int_seq_len){
        mask= (1<<8) - 1;
        n = uint8_actg(int_seq[p]) ;
        n4 = (n & mask);
        mask = mask <<8;
        n3 = (n & mask)>>8;
        mask = mask <<8;
        n2 = (n & mask)>>16;
        mask = mask <<8;
        n1 = (n & mask)>>24;

        actg_seq[p<<2] = n1;
        actg_seq[(p<<2) + 1] = n2;
        actg_seq[(p<<2) + 2] = n3;
        actg_seq[(p<<2) + 3] = n4;
        p = p + 1;
    }
}

typedef struct {
    int p; // suffix length; at least 8
    kc_c4_t **h; // 1<<p hash tables
} kc_c4x_t;

static kc_c4x_t *c4x_init(int p)
{
    int i;
    kc_c4x_t *h;
    CALLOC(h, 1);
    MALLOC(h->h, 1<<p);
    h->p = p;
    for (i = 0; i < 1<<p; ++i)
        h->h[i] = kc_c4_init();
    return h;
}

typedef struct {
    int n, m;
    uint64_t *a;
} buf_c4_t;


static inline void c4x_insert_buf(buf_c4_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
    int pre = y & ((1<<p) - 1);
    buf_c4_t *b = &buf[pre];
    if (b->n == b->m) {
        b->m = b->m < 8? 8 : b->m + (b->m>>1);
        REALLOC(b->a, b->m);
    }
    b->a[b->n++] = y;
}

static void count_seq_buf(buf_c4_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
    int i, l;
    uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
    for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) { // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
            if (++l >= k) { // we find a k-mer
                uint64_t y = x[0] < x[1]? x[0] : x[1];
                c4x_insert_buf(buf, p, hash64(y, mask));
            }
        } else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
    }
}

typedef struct { // global data structure for kt_pipeline()
    int k, block_len, n_thread;
    kseq_t *ks;
    kc_c4x_t *h;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
    pldat_t *p;
    int n, m, sum_len, nk;
    int *len;
    char **seq;
    buf_c4_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
    stepdat_t *s = (stepdat_t*)data;
    buf_c4_t *b = &s->buf[i];
    kc_c4_t *h = s->p->h->h[i];

    int j, p = s->p->h->p;
    for (j = 0; j < b->n; ++j) {
        khint_t k;
        int absent;
        k = kc_c4_put(h, b->a[j]>>p<<KC_BITS, &absent);
        if ((kh_key(h, k)&KC_MAX) < KC_MAX) ++kh_key(h, k);
    }
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    pldat_t *p = (pldat_t*)data;
    if (step == 0) { // step 1: read a block of sequences
        int ret;
        stepdat_t *s;
        CALLOC(s, 1);
        s->p = p;
        while ((ret = kseq_read(p->ks)) >= 0) {
            int l = p->ks->seq.l;
            if (l < p->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
            MALLOC(s->seq[s->n], l);
            memcpy(s->seq[s->n], p->ks->seq.s, l);
            s->len[s->n++] = l;
            s->sum_len += l;
            s->nk += l - p->k + 1;
            if (s->sum_len >= p->block_len)
                break;
        }
        if (s->sum_len == 0) free(s);
        else return s;
    } else if (step == 1) { // step 2: extract k-mers
        stepdat_t *s = (stepdat_t*)in;
        int i, n = 1<<p->h->p, m;
        CALLOC(s->buf, n);
        m = (int)(s->nk * 1.2 / n) + 1;
        for (i = 0; i < n; ++i) {
            s->buf[i].m = m;
            MALLOC(s->buf[i].a, m);
        }
        for (i = 0; i < s->n; ++i) {
            count_seq_buf(s->buf, p->k, p->h->p, s->len[i], s->seq[i]);
            free(s->seq[i]);
        }
        free(s->seq); free(s->len);
        return s;
    } else if (step == 2) { // step 3: insert k-mers to hash table
        stepdat_t *s = (stepdat_t*)in;
        int i, n = 1<<p->h->p;
        kt_for(p->n_thread, worker_for, s, n);
        for (i = 0; i < n; ++i) free(s->buf[i].a);
        free(s->buf); free(s);
    }
    return 0;
}



static kc_c4x_t *count_file(const char *fn, int k, int p, int block_size, int n_thread)
{
    pldat_t pl;
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.n_thread = n_thread;
    pl.h = c4x_init(p);
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}

static kc_c4x_t *count_file2(const char *fn, kc_c4x_t *hh, int k, int p, int block_size, int n_thread)
{
    pldat_t pl;
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.n_thread = n_thread;
    //pl.h = (kc_c4x_t *)hh;
    pl.h = hh;
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}

static kc_c4x_t *count_seq(const char *fn, int k, int p, int block_size, int n_thread)
{
    pldat_t pl;
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.n_thread = n_thread;
    pl.h = c4x_init(p);
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}




int kmer_cov(uint64_t kmer, uint64_t mask, kc_c4x_t *h){

    int j, x, cov=0, a_key;
    //unsigned cov=0;
    uint64_t hash_key = hash64(kmer, mask);
    j = hash_key & ((1<<KC_BITS) - 1);
    if(kh_size(h->h[j]) < 1){return 0;}

    //debug_print("kmer_cov here A\n");

    hash_key = hash_key >> KC_BITS<< KC_BITS;

    //debug_print("kmer_cov here B\n");

    x = kc_c4_get(h->h[j], hash_key);

    //debug_print("kmer_cov here C, j=%d\t",j);
    //debug_print("x=%d\n",x);

    //cov = 1004797953 & 16383;
    //cov = 1;
    //debug_print("kmer_cov here E, cov=%d\n", cov);

    //if(x == 0){return 0;} //To avoid of segment fault
    // if(x == 0){if(kh_size(h->h[j]) < 1){return 0;}} //To avoid of segment fault
    if(kh_exist(h->h[j], x)){
        //fprintf(stdout, "kmer_cov here D\n");
        a_key = kh_key(h->h[j], x); 
        //fprintf(stdout, "kmer_cov here E, a_key=%d\n", a_key);
        //fprintf(stdout, "kmer_cov here E, KC_MAX=%d\n", KC_MAX);
        cov = a_key & KC_MAX;

    }else{
        //fprintf(stdout, "kmer_cov here E\n");
        return 0;
    }
   
    // fprintf(stdout, "kmer_cov here F\n");
    // fprintf(stdout, "kmer_cov here G, cov=%d\n", cov);
    return cov;
}


static inline void add_kmer(uint64_t y, uint64_t mask, kc_c4x_t *h) //add k-mer to hash set
{
    int p = h->p;
    uint64_t y_hash = hash64(y, mask);
    int pre = y_hash & ((1<<p) - 1);
    khint_t k;
//        fprintf(stdout, "work for j %d\t", i);
    int absent;
    k = kc_c4_put(h->h[pre], y_hash>>p<<KC_BITS, &absent);
    ++kh_key(h->h[pre], k);
    //if ((kh_key(h[pre], k)&KC_MAX) < KC_MAX) ++kh_key(h[pre], k);

}

int insert_kms(uint64_t *kms, uint64_t y, int km_num) // insert a k-mer $y to a %kms array, Lifu Song
{
    int i;
    for(i=0; i< km_num; i++){
       if(y == *(kms + i)){
            return km_num;
       }
    }
    *(kms + km_num) = y;
    return km_num + 1;
}



int kms_max_cov(uint64_t *kms, int km_num, uint64_t mask, kc_c4x_t *h){
    int i, max_cov=0;
    //fprintf(stdout, "func kms_max_cov\n");
    for(i=0; i<km_num; i++){
            //printf("");
            //fprintf(stdout, "uint64: %"PRIu64, *(kms+i));
            //fprintf(stdout, "\n");
            int cov = kmer_cov(*(kms+i), mask, h);
            if(cov > max_cov){max_cov=cov;}
    }
    //printf("km_num: %d\n", km_num);
    //printf("max cov: %d\n", max_cov);
    return max_cov;
}

int kms_min_cov(uint64_t *kms, int km_num, uint64_t mask, kc_c4x_t *h){
    int i, min_cov=10000;
    //fprintf(stdout, "func kms_max_cov\n");
    for(i=0; i<km_num; i++){
            //printf("");
            //fprintf(stdout, "uint64: %"PRIu64, *(kms+i));
            //fprintf(stdout, "\n");
            int cov = kmer_cov(*(kms+i), mask, h);
            //fprintf(stderr, "cov: %d\n", cov);
            if(cov < min_cov){min_cov=cov;}
    }
    //printf("km_num: %d\n", km_num);
    //printf("max cov: %d\n", max_cov);
    return min_cov;
}

int kms_ext_num(uint64_t *kms, int km_num, uint64_t mask, kc_c4x_t *h){
    int i, km_ex_num=0;
    //fprintf(stdout, "func kms_max_cov\n");
    for(i=0; i<km_num; i++){
            //printf("");
            //fprintf(stdout, "uint64: %"PRIu64, *(kms+i));
            //fprintf(stdout, "\n");
            int cov = kmer_cov(*(kms+i), mask, h);
            //fprintf(stderr, "cov: %d\n", cov);
            if(cov > 0){km_ex_num=km_ex_num + 1;}
    }
    //printf("km_num: %d\n", km_num);
    //printf("max cov: %d\n", max_cov);
    return km_ex_num;
}

int kms_avg_cov(uint64_t *kms, int km_num, uint64_t mask, kc_c4x_t *h){
    int i, avg_cov=0, k_num=1;
    //fprintf(stdout, "func kms_max_cov\n");
    for(i=0; i<km_num; i++){
            //printf("");
            //fprintf(stdout, "uint64: %"PRIu64, *(kms+i));
            //fprintf(stdout, "\n");
            int cov = kmer_cov(*(kms+i), mask, h);
            if(cov > 0){avg_cov = cov + avg_cov; k_num = k_num + 1;}
    }
    //printf("km_num: %d\n", km_num);
    //printf("max cov: %d\n", max_cov);
    if(k_num>0){
        avg_cov = (int) avg_cov/k_num;
    }else{
        avg_cov=0;
    }

    return avg_cov;
}


int seq_kmers(uint64_t *kms, int k, int len, const char *seq) // insert k-mers in $seq to $kms
{
    int i, l, km_num=0;
    uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
    for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) { // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            if (++l >= k) { // we find a k-mer
				*(kms + km_num) = x[0];
				km_num = km_num + 1;
            }
        } else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
    }
    return km_num;
}

typedef struct {
    uint64_t kmer_s;
    uint64_t kmer_t;
	int pos_s;
    int pos_t;
    char *name[100];
}var_location;

typedef struct { // data structure for file evaluation
    kc_c4x_t *h;
    kc_c4x_t *hr;
    kc_c4x_t *hr_pos;
    float error_rate;

    uint64_t *kms;
    int k;
    int p;
    var_location *var_locs;
} evaluation_t; 

typedef struct {
    struct path_node* pre_node;
    uint64_t kmer;
    unsigned int cov;
	int pos;
}path_node;


int path_node_num(path_node *term_node){
    // fprintf(stderr, "##### path_node_num() \n");
    path_node *p_node;
    int node_num=1;
    p_node = term_node;

    // fprintf(stderr, "hello pathnode while \n");

    while(p_node->pre_node){
        // fprintf(stderr, "hello pathnode while \n");
        node_num++;
        p_node = p_node->pre_node;
    }
    // fprintf(stderr, "##### path_node_num: %d \n", node_num);
    return node_num;
}


unsigned char* nodes_path(path_node *term_node, unsigned char *path_seq, int k){
    path_node *p_node = term_node;
    int path_len = path_node_num(p_node);

    // fprintf(stdout, "\n\nnodes_path function start here\n");
    // fprintf(stdout, "path node num: %d\n", path_len);

    int passed_nodes = 0;   
    // print_uint64_kmer(p_node->kmer, km_len);
    //p_node = p_node->pre_node; //Skipping the terminal kmer 
    while(p_node->pre_node && passed_nodes < k - 1){
        p_node = p_node->pre_node;
        // print_uint64_kmer(p_node->kmer, km_len);
        //*(path_seq + path_seq_p) = (unsigned char)nt4_seq_table[p_node->kmer & 3ULL]; 
        passed_nodes++;
    }


    // fprintf(stdout, "\n\nnnode after moving:\n");
    // print_uint64_kmer(p_node->kmer, km_len);

    int path_seq_p = path_len-passed_nodes-2;   
    //p_node = p_node->pre_node; //Skipping the terminal kmer 
    while(p_node->pre_node && path_seq_p >= 0){
        p_node = p_node->pre_node;
        // print_uint64_kmer(p_node->kmer, km_len);
        // fprintf(stdout, "nodes path function: %d\n", path_seq_p - 1);
        *(path_seq + path_seq_p - 1) = (unsigned char)nt4_seq_table[p_node->kmer & 3ULL]; 
        // fprintf(stdout, "nodes path function: %c\n", (unsigned char)nt4_seq_table[p_node->kmer & 3ULL]);
        path_seq_p--;
    }
    return path_len;
}

// Working on 20230131 Lifu Song
uint64_t* nodes_to_kms(path_node *term_node, int k){
    // fprintf(stderr, "##### nodes_to_kms() \n");
    path_node *p_node = term_node;
    uint64_t *kms;

    int path_len = path_node_num(p_node);
    // fprintf(stderr, "##### nodes_to_kms() path_len: %d\n", path_len);

    MALLOC(kms, path_len);

    int passed_nodes = 0;   
    *(kms + path_len -1) = p_node->kmer; 
    while(p_node->pre_node ){
        p_node = p_node->pre_node;
        passed_nodes++;
        *(kms + path_len - passed_nodes -1) = p_node->kmer; 
        
    }
    return kms;
}


int nodes_path_cov( evaluation_t *eva, path_node *term_node){
    // fprintf(stderr, "nodes_path_cov() function ...... \n");
    kc_c4x_t *h = eva->h;
    path_node *p_node = term_node;
    int k = eva->k;
    float error_rate = eva->error_rate;
    uint64_t mask = (1ULL<<k*2) - 1; //


    int path_len = path_node_num(p_node);
    int path_cov = kmer_cov(min_hash_key(p_node->kmer, k),mask,h);
    int passed_nodes = 0, p_node_cov = 0;   

    
    while(p_node->pre_node && passed_nodes < k - 1){
        p_node = p_node->pre_node;
        p_node_cov = kmer_cov(min_hash_key(p_node->kmer, k),mask,h);
        if(p_node_cov < path_cov){
            path_cov = p_node_cov;
        }
        passed_nodes++;
    }
    return path_cov;
}


float nodes_path_p_value( evaluation_t *eva, path_node *term_node){
    
    kc_c4x_t *h = eva->h;
    path_node *p_node = term_node;
    int k = eva->k;
    float error_rate = eva->error_rate;
    uint64_t mask = (1ULL<<k*2) - 1; //

    //fprintf(stderr, "dfasdfasdf\t");

    int path_len = path_node_num(p_node);
    int path_cov = kmer_cov(min_hash_key(p_node->kmer, k),mask,h);
    int passed_nodes = 0, p_node_cov = 0;   

    

    while(p_node->pre_node && passed_nodes < k - 1){
        p_node = p_node->pre_node;
        p_node_cov = kmer_cov(min_hash_key(p_node->kmer, k),mask,h);
        if(p_node_cov < path_cov){
            path_cov = p_node_cov;
        }
        passed_nodes++;
    }

    float p_kmer = pow(1-error_rate,2*k) * pow(error_rate/3.0, path_len-k-1);
    float x, y;
    // x = pow(1-error_rate,2*k);
    // y = pow(error_rate/3.0, path_len-k-1);

    // fprintf(stdout, "\tx:%f\t", x);
    // fprintf(stdout, "\ty:%f\t", y);
    // fprintf(stdout, "\tcov:%d\t", path_cov);

    float p_value = pow(p_kmer,path_cov);
    //fprint(stderr, "dfasdfasdf\n");
    return p_value;
}




int kms_to_seq(unsigned char *aseq, uint64_t *kms, int ss, int tt){
    
    int i=0, l = tt-ss+1;
    //CALLOC(aseq, l);
    // fprintf(stdout, "ss: %d\t", ss);
    // fprintf(stdout, "tt: %d\n", tt);

    while(i<l){
        //print_uint64_kmer(kms[ss + i], 21);
        int n = kms[ss + i] % 4;
        *(aseq + i)  = nt4_seq_table[n];
        // fprintf(stdout, "nt4 seq table: %d\n", (int) kms[i] & 3ULL);
        //fprintf(stdout, "%d\t", n);
        //fprintf(stdout, "%c\n", nt4_seq_table[n]);
        i = i + 1;
    }
    aseq[l] = '\0';
    return aseq;
}

path_node all_nodes[Max_Path_Num * 2000];
path_node *term_nodes[Max_Path_Num+1]; 
path_node *good_term_nodes[Max_Path_Num+1]; 

// 
int slim_path(evaluation_t *eva, var_location *a_var, var_location *new_var, uint64_t *path_kms, int path_nodes_num){
    // fprintf(stderr, "##### slim_path() \n");
    

    //int slim_path_len = path_len;
    //fprintf(stderr, "##### slim_path_len: %d \n", path_len);
    // MALLOC(slim_kms, slim_path_len); 

    //var_location new_var;
    new_var->pos_s = 0;
    new_var->pos_t = path_nodes_num - 1;
 
    int pad_l = 0, pad_r = 0;

    //Slim right
    while(eva->kms[a_var->pos_t - pad_r] == path_kms[path_nodes_num - 1 - pad_r] &&  path_nodes_num - pad_l - pad_r > eva->k && a_var->pos_t - a_var->pos_s - pad_l - pad_r > eva->k){
    //while(eva->kms[a_var->pos_t - pad_r] == path_kms[path_nodes_num - 1 - pad_r]){
        pad_r++;
        // fprintf(stdout, "##### pad_r: %d \n", pad_r);
    }
    
    //fprintf(stderr, "##### Hello \n");
    //Slim left
    while(eva->kms[a_var->pos_s + pad_l] == path_kms[pad_l]                      && path_nodes_num - pad_l - pad_r > eva->k && a_var->pos_t - a_var->pos_s - pad_l - pad_r > eva->k ){
    //while(  eva->kms[a_var->pos_s + pad_l] == path_kms[pad_l] ){
        pad_l++;
        // fprintf(stdout, "##### pad_l: %d \n", pad_l);
    }
    //if(pad_l > 0){pad_l--;}
    
    if(pad_r > 0){pad_r--;}
    if(pad_l > 0){pad_l--;}

    
    new_var->pos_s = new_var->pos_s + pad_l;
    new_var->pos_t = new_var->pos_t - pad_r;

    a_var->pos_s = a_var->pos_s + pad_l;
    a_var->pos_t = a_var->pos_t - pad_r;

    if(pad_r > 0 || pad_l > 0){
        return 1;
    }else{
        return 0;
    } 
}


// 
int output_path(evaluation_t *eva, int var_loc_p, int path_index){ 
    // fprintf(stdout, "\n\n\n##### output_path() \n");
    kc_c4x_t *h = eva->h;
    int k = eva->k;
    uint64_t mask = (1ULL<<k*2) - 1; 

    //Sliming path
    var_location ref_var;
    ref_var.pos_s = eva->var_locs[var_loc_p].pos_s;
    ref_var.pos_t = eva->var_locs[var_loc_p].pos_t;

    // int ref_pos_s = eva->var_locs[var_loc_p].pos_s;
    // int ref_pos_t = eva->var_locs[var_loc_p].pos_t;

    path_node *p_node = good_term_nodes[path_index];
    int path_nodes_num = path_node_num(good_term_nodes[path_index]);
    
    uint64_t *path_kms = nodes_to_kms(p_node, k);
    // uint64_t *slim_kms;
    // fprintf(stdout, "##### path_nodes_num: %d\n", path_nodes_num);

    var_location seq_var;
    //paded = if pad happend or not. 
    // fprintf(stdout, "before slim, path node num: %d\n", path_nodes_num);
    int paded = slim_path(eva, &ref_var, &seq_var, path_kms, path_nodes_num);

    // fprintf(stdout, "##### slim_path function finished\n");

    int start_pos = ref_var.pos_s; 
    int term_pos = ref_var.pos_t; 

    // fprintf(stdout, "##### ref kmers\n");
    // print_uint64_kmer(eva->var_locs[var_loc_p].kmer_s, k); 
    // fprintf(stdout, "\t");
    // print_uint64_kmer(eva->var_locs[var_loc_p].kmer_t, k); 
    // fprintf(stdout, "\n");

    // fprintf(stdout, "##### var kmers\n");
    // print_uint64_kmer(ref_var.pos_s, k); 
    // print_uint64_kmer(ref_var.pos_t, k); 

    // fprintf(stdout, "##### ref_var.pos_s: %d \t", ref_var.pos_s);
    // fprintf(stdout, " ref_var.pos_t: %d \t", ref_var.pos_t);
    // fprintf(stdout, " new_var start_pos: %d \t", seq_var.pos_s);
    // fprintf(stdout, " new_var term_pos: %d \n", seq_var.pos_t);
    
    // fprintf(stdout, "#Outputting the varation analysis results\n");

    int path_cov = nodes_path_cov(eva, good_term_nodes[path_index]);

    // fprintf(stdout, "##### nodes_path_cov() function finished\n");

    unsigned char *ref_seq;
    int ref_seq_len = ref_var.pos_t - ref_var.pos_s - k + 2;
    CALLOC(ref_seq, ref_seq_len + 10);
    //kms_to_seq(ref_seq, eva->kms, ref_var.pos_s, ref_var.pos_t - k + 1 );

    unsigned char *path_seq;
    int slim_path_len = seq_var.pos_t - seq_var.pos_s - k + 2;
    CALLOC(path_seq, slim_path_len + 500); 
    //kms_to_seq(path_seq, path_kms, seq_var.pos_s, seq_var.pos_t - k + 1 );
    
    // fprintf(stdout, "##### Here\n");
    //slim_path_len = slim_path_len - k - 1;
    // fprintf(stdout, "slim_path_len: %d\t", slim_path_len);
    // fprintf(stdout, "ref_seq_len: %d\n", ref_seq_len);

    fprintf(stdout, "%s\t",eva->var_locs[var_loc_p].name);
    if(slim_path_len >= ref_seq_len){
        kms_to_seq(ref_seq, eva->kms, ref_var.pos_s+1, ref_var.pos_t - k  );
        kms_to_seq(path_seq, path_kms, seq_var.pos_s+1, seq_var.pos_t - k  );
        fprintf(stdout, "%d\t", ref_var.pos_s + k + 1); 
        fprintf(stdout, ".\t");
        fprintf(stdout, "%s\t", ref_seq);
        fprintf(stdout, "%s\t", path_seq);
        fprintf(stdout, "%d\t", path_cov);
        if(slim_path_len > ref_seq_len){fprintf(stdout, "INS\t"); }else{fprintf(stdout, "SUB\t"); }

    }
    if(slim_path_len < ref_seq_len){ //Deletion
        // ref_var.pos_s = ref_var.pos_s -1;
        // seq_var.pos_s = seq_var.pos_s -1;
        // ref_seq_len++; slim_path_len++;
        if(slim_path_len < 0){
            fprintf(stderr, "Slim Path less than 0 detected\n");
            kms_to_seq(ref_seq, eva->kms, ref_var.pos_s + slim_path_len, ref_var.pos_t - k + 1 );
        }else{
            kms_to_seq(ref_seq, eva->kms, ref_var.pos_s , ref_var.pos_t - k + 1 );
            kms_to_seq(path_seq, path_kms, seq_var.pos_s, seq_var.pos_t - k + 1 );
        }
        //kms_to_seq(ref_seq, eva->kms, ref_var.pos_s, ref_var.pos_t - k + 1 );
        fprintf(stdout, "%d\t", ref_var.pos_s + k );
        fprintf(stdout, ".\t");
        fprintf(stdout, "%s\t", ref_seq);
        fprintf(stdout, "%s", path_seq); 
        fprintf(stdout, "\t%d\t", path_cov);
        fprintf(stdout, "DEL\t");
    }
    //fprintf(stdout, "\n");
    
    // fprintf(stdout, "##### ref_var.pos_s: %d \t", ref_var.pos_s);
    // fprintf(stdout, " ref_var.pos_t: %d \t", ref_var.pos_t);
    // fprintf(stdout, " new_var start_pos: %d \t", seq_var.pos_s);
    // fprintf(stdout, " new_var term_pos: %d \n\n\n\n\n", seq_var.pos_t);

}

// Greedy path search function
int var_path_search_ref(evaluation_t *eva, int var_loc_p, int max_path_len, int min_cov){ 
    // fprintf(stderr, "##### var_path_search_ref() \n");
	kc_c4x_t *h = eva->h;
    int k = eva->k;
    uint64_t mask = (1ULL<<k*2) - 1; //
    var_location a_var = eva->var_locs[var_loc_p];
    uint64_t kmer_s = a_var.kmer_s;
    uint64_t kmer_t = a_var.kmer_t;
    int start_pos = a_var.pos_s; 
    int term_pos = a_var.pos_t; 
    
    /* debuging scripts

        debug_print(stdout, "\n\n\n\n\n#####\nVaration location: %d\n", var_loc_p);

        debug_print(stdout, "\n%d\t", start_pos);
        debug_print(stdout, "%d\n", term_pos);

        debug_print(stdout, "\nStart k-mer:\n");
        print_uint64_kmer(kmer_s, k);
        fprintf(stdout, "\nStart k-mer cov: %d\n", kmer_cov(min_hash_key(kmer_s, k), mask,h));


        debug_print(stdout, "\nTerm k-mer:\n");
        print_uint64_kmer(kmer_t, k); 
        debug_print(stdout, "\nTerm k-mer cov: %d\n", kmer_cov(min_hash_key(kmer_t, k), mask,h));
    
    */
    
    all_nodes[0].kmer = kmer_s;
    all_nodes[0].cov = kmer_cov(min_hash_key(kmer_s, k),mask,h);

    int a_cov = 1;
    if(a_cov < min_cov ){a_cov = min_cov ;}

    term_nodes[0] = &all_nodes[0]; 

    int fresh_terms = 1, all_node_num = 1, p = 0;

    path_node *p_node;
    path_node next_nodes[4];
    next_nodes[0].cov = 0; next_nodes[1].cov = 0; next_nodes[2].cov = 0; next_nodes[3].cov = 0;

    path_node *pre_term_nodes[Max_Path_Num+1]; 
    int good_term_node_num = 0;

    while(p <= max_path_len && fresh_terms > 0 && fresh_terms <= Max_Path_Num){
        // fprintf(stdout, "Hello A\n");
        int i;
        //Copying old term nodes pointer
        for(i=0; i < fresh_terms; i++){pre_term_nodes[i] = term_nodes[i]; }
        //new_term_node_num = 0;
        int new_fresh_terms = 0;
        int n_p = 0;
        while(n_p < fresh_terms && new_fresh_terms <=Max_Path_Num ){
            p_node = pre_term_nodes[n_p]; //Obtain term node one by one
            // fprintf(stdout, "Hello B\n");
            // else{
            //a_cov = p_node->cov/Max_Cov_Ratio;
            //if(a_cov < cut_cov ){a_cov = cut_cov;}
            uint64_t kp = (*p_node).kmer<<2;

            next_nodes[0].kmer = kp & mask;
            next_nodes[1].kmer = (kp + 1) & mask;
            next_nodes[2].kmer = (kp + 2) & mask;
            next_nodes[3].kmer = (kp + 3) & mask;
            next_nodes[0].cov = kmer_cov(min_hash_key(next_nodes[0].kmer, k),mask,h);
            next_nodes[1].cov = kmer_cov(min_hash_key(next_nodes[1].kmer, k),mask,h);
            next_nodes[2].cov = kmer_cov(min_hash_key(next_nodes[2].kmer, k),mask,h);
            next_nodes[3].cov = kmer_cov(min_hash_key(next_nodes[3].kmer, k),mask,h);

            // fprintf(stdout, "cov: %d\t", next_nodes[0].cov);
            // fprintf(stdout, "%d\t", next_nodes[1].cov);
            // fprintf(stdout, "%d\t", next_nodes[2].cov);
            // fprintf(stdout, "%d\n", next_nodes[3].cov);

            int i;
            for(i=0;i<4;i++){
                // if(next_nodes[i].cov >=a_cov && next_nodes[i].cov <= b_cov){
                if(next_nodes[i].cov >=a_cov){
                    //print_uint64_kmer(next_nodes[i].kmer, k);
                    //print_uint64_kmer(kmer_t, k);
                    if( next_nodes[i].kmer == kmer_t){ // if the k-mer is equal to kmer_t, then 
                      // skip this term node during path searching  
                        all_nodes[all_node_num].pre_node = p_node;
                        all_nodes[all_node_num].kmer = next_nodes[i].kmer;
                        all_nodes[all_node_num].cov = next_nodes[i].cov;
                        good_term_nodes[good_term_node_num] = &all_nodes[all_node_num];
                        good_term_node_num++; 
                    }else{
                        all_nodes[all_node_num].pre_node = p_node;
                        all_nodes[all_node_num].kmer = next_nodes[i].kmer;
                        all_nodes[all_node_num].cov = next_nodes[i].cov;
                        term_nodes[new_fresh_terms] = &all_nodes[all_node_num];
                        //new_term_node_num++;
                        new_fresh_terms++; 
                    }
                    all_node_num++;
                }  
            }
            n_p++;
        }
        fresh_terms = new_fresh_terms;
        p++;
    }
    if(good_term_node_num == 1){
        //good_term_nodes
        if(path_node_num(good_term_nodes[good_term_node_num-1]) < k ){
            good_term_node_num = 0;
        }
    }
    // try to search from term node
    return good_term_node_num; 
}

// Assembly-based variation analysis of sepecific location. The greedy path searching is performed by var_path_search_ref() function.
int var_analysis_ref(evaluation_t *eva, int var_loc_p, int max_path_len, int min_cov){ 
    // fprintf(stderr, "##### var_analysis_ref() \n");
    // debuging scripts
    uint64_t mask = (1ULL<<eva->k*2) - 1;
    debug_print("\n\n\n\n\n#####\nVaration location: %d\n", var_loc_p);

    debug_print("\n%d\t", eva->var_locs[var_loc_p].pos_s);
    debug_print("%d\n", eva->var_locs[var_loc_p].pos_t);

    debug_print("\nStart k-mer:\n");
    print_uint64_kmer(eva->var_locs[var_loc_p].kmer_s, eva->k);
    debug_print("\nStart k-mer cov: %d\n", kmer_cov(min_hash_key(eva->var_locs[var_loc_p].kmer_s, eva->k), mask, eva->h));


    debug_print("\nTerm k-mer:\n");
    print_uint64_kmer(eva->var_locs[var_loc_p].kmer_t, eva->k); 
    debug_print("\nTerm k-mer cov: %d\n", kmer_cov(min_hash_key(eva->var_locs[var_loc_p].kmer_t, eva->k), mask, eva->h));
    
    

    // greedy path search
    var_location ori_var;
    ori_var.kmer_s = eva->var_locs[var_loc_p].kmer_s;
    ori_var.kmer_t = eva->var_locs[var_loc_p].kmer_t;
    ori_var.pos_s = eva->var_locs[var_loc_p].pos_s;
    ori_var.pos_t = eva->var_locs[var_loc_p].pos_t;


    int good_term_node_num = var_path_search_ref(eva, var_loc_p, max_path_len, min_cov);
    int ext_times = 0 ;
    while( good_term_node_num < 1 && ext_times < 5){ // if no path obtained, extent the node position and search again.
        extent_var_loc(eva, eva->kms, &eva->var_locs[var_loc_p], min_cov, eva->k);   
        good_term_node_num = var_path_search_ref(eva, var_loc_p, max_path_len, min_cov);
        ext_times = ext_times + 1;
    }
    
    int best_termnode_index = best_term_node(eva, eva->var_locs[var_loc_p].pos_t - eva->var_locs[var_loc_p].pos_s, good_term_node_num);
    output_path(eva, var_loc_p, best_termnode_index);
    if(good_term_node_num > 1){
        fprintf(stdout, "\tMultiple-path");
    }
    fprintf(stdout, "\n");

    //fprintf(stderr, "##### Good termnode num: %d \n", good_term_node_num);
    // int i;
    // for(i=0; i < good_term_node_num; i++){ output_path(eva, var_loc_p, i); }
    // // Obtaining var information

    // if(good_term_node_num < 1){
    //     fprintf(stdout, "%s\t",eva->var_locs[var_loc_p].name);
    //     fprintf(stdout, "%d\t", ori_var.pos_s + eva->k + 1); // +1 since the start base is pos 0
    //     fprintf(stdout, "\t%d\tFailed to detect the variation\n", eva->var_locs[var_loc_p].pos_t);
    // }

    return good_term_node_num; 
}

int best_term_node(evaluation_t *eva, int ref_node_num, int good_term_node_num){

    int k = eva->k;
    int i = 0, best_path=0, short_path_var_len = 10000, short_path_cov = 1;
    int p_path_node_num, var_len=0, p_path_cov;

    for(i=0; i < good_term_node_num; i++){
        p_path_node_num = path_node_num(good_term_nodes[i]);
        var_len = abs(p_path_node_num -ref_node_num) ;
        p_path_cov = nodes_path_cov(eva, good_term_nodes[i]);
        if(var_len < short_path_var_len || p_path_node_num == short_path_var_len && p_path_cov > short_path_cov){
            best_path = i;
            short_path_var_len = var_len;
            short_path_cov = p_path_cov;
        }
    }
    return best_path;
}

int extent_var_loc(evaluation_t *eva, uint64_t *kms, var_location *var_loc, int min_cov, int dis){ 
    
    // fprintf(stdout, "\nExten_var_loc\n");
    int k = eva->k;
    uint64_t mask = (1ULL<<k*2) - 1;

    var_loc->pos_s = var_loc->pos_s - dis;
    var_loc->pos_t = var_loc->pos_t + dis;

    int cov_s = kmer_cov(min_hash_key(kms[var_loc->pos_s],k), mask, eva->h);
    int cov_t = kmer_cov(min_hash_key(kms[var_loc->pos_t],k), mask, eva->h);
    //int rep_km_num_s = kmer_cov(min_hash_key(kms[var_loc->pos_s],k), mask, eva->hr); 
    //int rep_km_num_t = kmer_cov(min_hash_key(kms[var_loc->pos_t],k), mask, eva->hr); 

    int itt = 0;
    while(cov_s <= min_cov && itt < dis){
        var_loc->pos_s = var_loc->pos_s - 1;
        cov_s = kmer_cov(min_hash_key(kms[var_loc->pos_s],k), mask, eva->h);
        itt = itt + 1;
    }

    itt = 0;
    while(cov_t <= min_cov && itt < dis){
        var_loc->pos_t = var_loc->pos_t + 1;
        cov_t = kmer_cov(min_hash_key(kms[var_loc->pos_t],k), mask, eva->h);
        itt = itt + 1;
    }

    var_loc->kmer_s = kms[var_loc->pos_s];
    var_loc->kmer_t = kms[var_loc->pos_t];
    
    int cov_rr = 0, cov_tt =0;
    cov_rr = kmer_cov(min_hash_key(var_loc->kmer_s,k), mask, eva->hr);
    cov_tt = kmer_cov(min_hash_key(var_loc->kmer_t,k), mask, eva->hr);

    if(cov_rr >1 || cov_tt > 1){
        extent_var_loc(eva, kms, var_loc, min_cov, dis); 
    }

    return 1;
} 


int optimize_var_loc(evaluation_t *eva, uint64_t *kms, var_location *var_loc, int min_cov){ 
    //fprintf(stdout, "#############optimize_var_loc() function \n");
    int k = eva->k;
    //int p = eva->p;
    uint64_t mask = (1ULL<<k*2) - 1;
    int cov_s = kmer_cov(min_hash_key(kms[var_loc->pos_s],k), mask, eva->h);
    int cov_t = kmer_cov(min_hash_key(kms[var_loc->pos_t],k), mask, eva->h);
    int rep_km_num = kmer_cov(min_hash_key(kms[var_loc->pos_s],k), mask, eva->hr); 

    // fprintf(stdout, "Optimizing start posA: %d \n", var_loc->pos_s);
    // fprintf(stdout, "Cov_s %d \n", cov_s);
    // fprintf(stdout, "Cov_t: %d \n", cov_t);

    int shift_time = 0;
    while(cov_s > min_cov && cov_t/cov_s > 2.8 || rep_km_num > 2 && shift_time <= k){
        cov_s = kmer_cov(min_hash_key(kms[var_loc->pos_s - 1],k), mask, eva->h);
        if(cov_s >= min_cov && cov_s > 0){
            var_loc->pos_s = var_loc->pos_s -1;
            var_loc->kmer_s = kms[var_loc->pos_s];
            
            rep_km_num = kmer_cov(min_hash_key(kms[var_loc->pos_s],k), mask, eva->hr); 
            // fprintf(stdout, "Optimizing start posB: %d \n", var_loc->pos_s);
            // fprintf(stdout, "Cov_s %d \n", cov_s);
            // fprintf(stdout, "Cov_t: %d \n", cov_t);
            shift_time = shift_time + 1;

        }else{cov_s = 0;}
        
    }
    return 1;
} 

// Output the reference seq k-mer coverages
static evaluation_t *analysis_ref_seq(evaluation_t *eva, const char *fn, int max_len, int min_cov)
{
    pldat_t pl;
    gzFile fp;
    int k = eva->k;
    int p = eva->p;

    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = eva->k;

    fprintf(stderr, "Analyzing the locations with variations ...... \n"); 
    uint64_t mask = (1ULL<<k*2) - 1;
    int ret;
    stepdat_t *s;
    CALLOC(s, 1);
    s->p = &pl;
    int i = 0;
    var_location *var_locs;
    
    while ((ret = kseq_read(pl.ks)) >= 0) {// Reading one seq for each loop
        //eva->strand_num++;
        
        i = i + 1;
        uint64_t *kms;  //all k-mers in array, all forward direction

        int l = pl.ks->seq.l, km_num = 0, miss_km_num = 0; 
        if (l < k) continue;

        MALLOC(kms, l); 
        km_num = seq_kmers(kms, k, l, pl.ks->seq.s);
        
        MALLOC(var_locs, (l-k + 1));

        int j = 0;
        int var_loc_num = 0;
        var_location one_var_loc;
        one_var_loc.kmer_s = 0;
        one_var_loc.pos_s = -1;
        //one_var_loc.name = pl.ks->name;
        pl.ks->name.s;

        //track type=wiggle_0 name="k-mer coverage" 
        //fixedStep chrom=NC_003450.3 start=1 step=1

        fprintf(stdout, "track type=wiggle_0 name=\"k-mer coverage\"\nfixedStep chrom=");
        //fprintf(stdout, "\"");
        fprintf(stdout, "%s start=1 step=1\n", pl.ks->name.s); 
        
        //uint64_t hs_value;
        int cov=0;
        int i=1;
        while(i < k){  //outputting the initial coverage values
            fprintf(stdout, "0.0\n"); 
            i++;
        }
        while(j < km_num){ // Output the k-mer coverage values
            cov = kmer_cov( min_hash_key(kms[j],k), mask, eva->h);
            fprintf(stdout, "%d.0\n", cov); 
            j++;
        } //var_locations obtained.  


  
    }

    
    // Cleaning 
    kseq_destroy(pl.ks);
    // gzclose(fp);
    // fprintf(stdout, "I am here, working .... \n");  
    
    return eva;
}


static kc_c4x_t *index_ref(const char *fn, int k, int p)
{
    pldat_t pl;
    gzFile fp;
    kc_c4x_t *h;
    pl.h = h;
    //c4x_init(p);

    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;

    //fprintf(stderr, "Analyzing the locations with variations ...... \n"); 
    //fprintf(stderr, "Analyzing the locations with variations ...... \n"); 
    uint64_t mask = (1ULL<<k*2) - 1;
    
    int ret;
    stepdat_t *s;
    CALLOC(s, 1);
    s->p = &pl;
    int sn = 0;
    
    while ((ret = kseq_read(pl.ks)) >= 0) {// Reading one seq for each loop
        //eva->strand_num++;
        sn = sn + 1;
        int i, l, len = pl.ks->seq.l;
        char *seq = pl.ks->seq.s;
        if (l < k) continue;
        // index k-mers
        kc_c4_t *hh ;
        uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
        for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
            int c = seq_nt4_table[(uint8_t)seq[i]];
            if (c < 4) { // not an "N" base
                x[0] = (x[0] << 2 | c) & mask;                  // forward strand
                //x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
                if (++l >= k) { // we find a k-mer
                    //uint64_t y = x[0] < x[1]? x[0] : x[1];
                    //uint64_t y = x[0];
                    khint_t k;
                    int absent;
                    uint64_t y = hash64(x[0], mask);
                    int pre = y & ((1<<p) - 1);
                    hh = &h[pre]; 
                    k = kc_c4_put(hh, x[0]>>p<<KC_BITS, &absent);
                    if ((kh_key(hh, k)&KC_MAX) < KC_MAX) kh_key(hh, k) = kh_key(hh, k) + i + 1;
                }
            } else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
        }
    }

    // Cleaning 
    kseq_destroy(pl.ks);
    gzclose(fp);
    
    return h;
}

static void worker_for_test(void *data, long i, int tid) // callback for kt_for()
{
    stepdat_t *s = (stepdat_t*)data;
    buf_c4_t *b = &s->buf[i];
    kc_c4_t *h = s->p->h->h[i];

    int j, p = s->p->h->p;
    for (j = 0; j < b->n; ++j) {
        khint_t k;
        int absent;
        k = kc_c4_put(h, b->a[j]>>p<<KC_BITS, &absent);
        if ((kh_key(h, k)&KC_MAX) < KC_MAX) ++kh_key(h, k);
    }
}


void print_uint64_kmer(uint64_t km, int km_len){ // Function for debugging. 
    //fprintf(stdout, "uint64: %"PRIu64, km);
   //fprintf(stdout, " rev: %"PRIu64, comp_rev(km, km_len));
   //fprintf(stdout, " minHash: %"PRIu64, min_hash_key(km,km_len));
    unsigned char *seq;
    MALLOC(seq, km_len+1);
    seq[km_len+1] = '\0';
    seq[km_len] = '\0';
    uint64_acgt(km, seq, km_len);
    //fprintf(stdout, "\tkm: %s\n", seq);
    debug_print("%s", seq);
    free(seq);
}

void print_kmer_seq(unsigned char *seq, int k){ //Debugging function
    int i=0;
    for(i=0;i<k;i++){
            printf("%c", *(seq + i));
    }
}

void print_cov(unsigned char* seq, int k, kc_c4x_t *h){ //Debugging function
    uint64_t mask = (1ULL<<k*2) - 1;
    int test_cov=0;
    test_cov = kmer_cov( actgkmer_hashkey(seq, k), mask, h);
    print_kmer_seq(seq, k);
    printf("\t");
    printf("cov: %d\n", test_cov);
}

int kc_c4x_t_size(kc_c4x_t *h, int p){
    int total_size = 0, i;
    for(i = 0; i < 1<<p; ++i) {
        total_size = total_size + kh_size(h->h[i]);
    }
    return total_size;
}

int kc_c4x_t_kmers(kc_c4x_t *h, int p, int cov){
    int j, total_size = 0;
    if(cov < 1){return kc_c4x_t_size(h,p);}
    for (j=0; j < 1<<p; ++j){
        kc_c4_t *g = h->h[j];
        khint_t kk;
        for (kk = kh_end(g); kk > 0; --kk)
            if (kh_exist(g, kk)){
                int c = kh_key(g, kk) & KC_MAX;
                //c = c < 255? c : 255;
                if (c > cov){ // delete kmers
                   total_size++;
                }
            }
    }
    return total_size;
}

void usage(int k, int n_thread, int min_cov, int insert_size, float error_rate){

        //fprintf(stderr, "\n=======================================================================================================\n");
        fprintf(stderr, "                                                                                               \n");
        fprintf(stderr, "A k-mer coverage analysis tool \n\(Version 20230301");
        fprintf(stderr, " Author: Lifu Song songlf@tib.cas.cn\) \n\n\n");

        fprintf(stderr, "This tool counts the k-mers in the NGS files and compare the results with a reference genome. \n");
        fprintf(stderr, "The coverages of each k-mers in the reference genome will be outputted in *.wig format, which can be directly loaded into IGV for visulization.\n");
        //fprintf(stderr, "The coverages of each k-mers in the reference genome will be outputted in *.wig format, \n");
        
        

        fprintf(stderr, "                                                                                 ");
        //fprintf(stderr, "\n=======================================================================================================\n\n");
        fprintf(stderr, "\n\n");
        fprintf(stderr, "Usage: \n\nkcov [options] Reference NGS_files > kmer-coverage.wig \n");
        fprintf(stderr, "                       [Supporting formats: *fq, *fa, *fq.gz, *fa.gz]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
        //fprintf(stderr, "  -c INT     minimal k-mer coverage for variant calling [%d]\n", min_cov);
        //fprintf(stderr, "  -l INT     maximal assembly length [%d]\n", insert_size);
        fprintf(stderr, "  -t INT     number of threads for k-mer counting (default [%d])\n", n_thread);
        
        fprintf(stderr, "\nSample commonds:\n ");
        fprintf(stderr, "\nkcov  ref_genome.fa reads.fq  > kmer-coverage.wig\n ");
        fprintf(stderr, "\nkcov  ref_genome.fa reads1.fq reads2.fq > kmer-coverage.wig\n ");
        fprintf(stderr, "\nkcov  ref_genome.fa reads1.fq reads2.fq reads3.fq > kmer-coverage.wig\n ");
        //fprintf(stderr, "  -e INT     percentages of sequencing error rate (for P-value calculation, default [%g])\n", error_rate);
        fprintf(stderr, "\n");
        return 0;
}

int main(int argc, char *argv[])
{
    int i, c, k = 31, p = KC_BITS, block_size = 10000000, n_thread = 3, min_cov = 8, insert_size = 500;
    float error_rate = 0.025f;

    ketopt_t o = KETOPT_INIT;
    while ((c = ketopt(&o, argc, argv, 1, "k:t:c:l:e", 0)) >= 0) {
        if (c == 'k') k = atoi(o.arg);
        //else if (c == 'p') p = atoi(o.arg);
        else if (c == 't') n_thread = atoi(o.arg);
        else if (c == 'c') min_cov = atoi(o.arg);
        else if (c == 'l') insert_size = atoi(o.arg); 
        else if (c == 'e') error_rate = atof(o.arg);
         
    }

    if (argc - o.ind < 2) {
        usage(k, n_thread, min_cov, insert_size, error_rate);
        return 1;
    }


    fprintf(stderr, "k-mer size: %d\n", k);
    fprintf(stderr, "Counting k-mers of NGS file 1  ......\n");

    kc_c4x_t *h, *hr, *hr_pos;
    h = count_file(argv[o.ind + 1], k, p, block_size, n_thread);
    

    int f_num = argc - o.ind, c_f_n = 2;
    while(c_f_n < f_num ) {
        fprintf(stderr, "Counting NGS file %d ......\n", c_f_n );
        h = count_file2(argv[o.ind + c_f_n ], h, k, p, block_size, n_thread);  
        c_f_n = c_f_n + 1;
    }

    //fprintf(stderr, "Analyzing k-mers of reference sequence ......\n");
    //hr = count_file(argv[o.ind + 1], k, p, block_size, n_thread);

    fprintf(stderr, "Analyzing k-mers of reference sequence ......\n");
    hr = count_file(argv[o.ind], k, p, block_size, n_thread);

    evaluation_t eva;
    eva.h = h;
    eva.hr =hr;
    eva.hr_pos =hr_pos;
    eva.k=k;
    eva.p=p;
    eva.error_rate = error_rate;


    
    analysis_ref_seq(&eva, argv[o.ind], insert_size, min_cov);


//Cleaning memories....
    for(i = 0; i < 1<<p; ++i) {
        kc_c4_destroy(h->h[i]);
        //kc_c4_destroy(eva.h->h[i]);
    }

    free(h->h); free(h);
    // free(hr->h); free(hr);
    //free(&eva); 
    return 0;
}

