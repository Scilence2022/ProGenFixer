#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h> 
#include <math.h>
#include <unistd.h>

#include "ketopt.h" // command-line argument parser
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop


#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

// Forward declaration to avoid implicit declaration warning for usage()
void usage(int k, int n_thread, int min_cov, int assem_min_cov, int insert_size, float error_rate, int max_assem_cov);
#include "khashl.h" // hash table
#define KC_BITS 12
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
    tmp = key ^ key >> 28; 	key = key ^ tmp >> 28; 	tmp = key ^ tmp >> 28; 	key = key ^ tmp >> 28;
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
    uint64_t x[2], mask = (1ULL<<k*2) - 1;
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
    char name[100];  // Change from pointer array to character array
} var_location;

// Add new struct for variations
typedef struct {
    char chrom[100];
    int pos;
    char ref[1000];
    char alt[1000];
    char type[10];
} variation_t;

// Add to evaluation_t struct
typedef struct { // data structure for file evaluation
    kc_c4x_t *h;
    kc_c4x_t *hr;
    kc_c4x_t *hr_pos;
    float error_rate;
    int fix_enabled;  // New flag for fix option
    variation_t *variations;
    int var_count;
    int var_capacity;
    int max_assem_cov; // Added max assembly coverage

    uint64_t *kms;
    int k;
    int p;
    var_location *var_locs;
    FILE *vcf_out;
    int iteration;
    char *output_base;  // Add this field to store the output base name
} evaluation_t;

// Forward declarations to avoid implicit declaration warnings for functions used before definition
int extent_var_loc(evaluation_t *eva, uint64_t *kms, var_location *var_loc, int assem_min_cov, int dis);
int best_term_node(evaluation_t *eva, int ref_node_num, int good_term_node_num);

typedef struct path_node {
    struct path_node* pre_node;
    uint64_t kmer;
    unsigned int cov;
    int pos;
} path_node;




typedef struct {
	uint64_t *c; // Changed from fixed array to pointer
} buf_cnt_t;

typedef struct {
	const kc_c4x_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_hist(void *data, long i, int tid) // callback for kt_for()
{
	hist_aux_t *a = (hist_aux_t*)data;
	uint64_t *cnt = a->cnt[tid].c; // Access pointer
	kc_c4_t *g = a->h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k)) {
			int c = kh_key(g, k) & KC_MAX;
			// Use c directly as the index, capped by KC_MAX
			if (c <= KC_MAX) { // Check added just in case, though c should always be <= KC_MAX
				++cnt[c];
			}
		}
}

static void print_hist(const kc_c4x_t *h, int n_thread, const char *output_base)
{
	hist_aux_t a;
	uint64_t *cnt; // Changed from fixed array to pointer
	int i, j;
	size_t hist_size = (size_t)KC_MAX + 1; // Calculate histogram size based on KC_MAX

	a.h = h;
	CALLOC(a.cnt, n_thread); // Allocate array of structs

	// Allocate count arrays within each struct
	for (j = 0; j < n_thread; ++j) {
		CALLOC(a.cnt[j].c, hist_size); // Allocate based on hist_size
		if (a.cnt[j].c == NULL) {
			fprintf(stderr, "Error: Memory allocation failed for histogram counts (thread %d)\n", j);
			// Cleanup partially allocated memory
			while (--j >= 0) free(a.cnt[j].c);
			free(a.cnt);
			return; // Exit if allocation fails
		}
	}

	kt_for(n_thread, worker_hist, &a, 1<<h->p);

	// Allocate the final aggregated count array
	CALLOC(cnt, hist_size);
	if (cnt == NULL) {
		fprintf(stderr, "Error: Memory allocation failed for aggregated histogram counts\n");
		// Cleanup worker thread allocations
		for (j = 0; j < n_thread; ++j) free(a.cnt[j].c);
		free(a.cnt);
		return; // Exit if allocation fails
	}

	// Aggregate counts from threads
	for (j = 0; j < n_thread; ++j) {
		for (i = 0; i < hist_size; ++i) { // Iterate up to hist_size
			cnt[i] += a.cnt[j].c[i];
		}
		free(a.cnt[j].c); // Free individual thread count arrays after aggregation
	}
	free(a.cnt); // Free the array of structs

	// Create output filename using the output_base
	char histo_fn[256];
	snprintf(histo_fn, sizeof(histo_fn), "%s.histo", output_base);
	
	// Open the output file
	FILE *fp = fopen(histo_fn, "w");
	if (!fp) {
		fprintf(stderr, "Error: Could not open histogram file %s\n", histo_fn);
		return;
	}
	
	// Write histogram data to file
	for (i = 1; i < hist_size; ++i) // Iterate up to hist_size
		fprintf(fp, "%d\t%ld\n", i, (long)cnt[i]);
	
	fclose(fp);
	fprintf(stderr, "K-mer histogram written to %s\n", histo_fn);

	// Free the aggregated count array
	free(cnt);
}


// Add comparison function
int compare_variations(const void *a, const void *b) {
    const variation_t *va = (const variation_t *)a;
    const variation_t *vb = (const variation_t *)b;
    return vb->pos - va->pos; // Sort descending by position
}

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


int nodes_path(path_node *term_node, unsigned char *path_seq, int k){
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
    uint64_t mask = (1ULL<<k*2) - 1; //


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
    float p_value = pow(p_kmer,path_cov);
    //fprint(stderr, "dfasdfasdf\n");
    return p_value;
}




unsigned char* kms_to_seq(unsigned char *aseq, uint64_t *kms, int start, int term){
    
    int i=0, l = term-start+1;
    //CALLOC(aseq, l);
    // fprintf(stdout, "ss: %d\t", ss);
    // fprintf(stdout, "tt: %d\n", tt);

    while(i<l){
        //print_uint64_kmer(kms[ss + i], 21);
        int n = kms[start + i] % 4;
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
    int k = eva->k;
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
    slim_path(eva, &ref_var, &seq_var, path_kms, path_nodes_num);

    // fprintf(stdout, "##### slim_path function finished\n");

    // Removed unused variables start_pos and term_pos

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

    fprintf(eva->vcf_out, "%s\t",eva->var_locs[var_loc_p].name);
    char info_field[200];
    char var_type[4];

    if(slim_path_len >= ref_seq_len){
        kms_to_seq(ref_seq, eva->kms, ref_var.pos_s+1, ref_var.pos_t - k  );
        kms_to_seq(path_seq, path_kms, seq_var.pos_s+1, seq_var.pos_t - k  );
        fprintf(eva->vcf_out, "%d\t", ref_var.pos_s + k + 1); 
        fprintf(eva->vcf_out, ".\t"); // ID
        fprintf(eva->vcf_out, "%s\t", ref_seq);
        fprintf(eva->vcf_out, "%s\t", path_seq);
        fprintf(eva->vcf_out, ".\t"); // QUAL
        fprintf(eva->vcf_out, "PASS\t"); // FILTER
        if(slim_path_len > ref_seq_len){strcpy(var_type, "INS"); }else{strcpy(var_type, "SUB"); }
        snprintf(info_field, sizeof(info_field), "KMER_COV=%d;VARTYPE=%s", path_cov, var_type);
        fprintf(eva->vcf_out, "%s", info_field);

    }
    if(slim_path_len < ref_seq_len){ //Deletion
        // ref_var.pos_s = ref_var.pos_s -1;
        // seq_var.pos_s = seq_var.pos_s -1;
        // ref_seq_len++; slim_path_len++;
        if(slim_path_len < 0){
            fprintf(stderr, "Slim Path less than 0 detected\n");
            kms_to_seq(ref_seq, eva->kms, ref_var.pos_s + slim_path_len + 1, ref_var.pos_t - k + 1 );
        }else{
            kms_to_seq(ref_seq, eva->kms, ref_var.pos_s+1, ref_var.pos_t - k  );
            kms_to_seq(path_seq, path_kms, seq_var.pos_s+1, seq_var.pos_t - k  );
        }
        //kms_to_seq(ref_seq, eva->kms, ref_var.pos_s, ref_var.pos_t - k + 1 );
        fprintf(eva->vcf_out, "%d\t", ref_var.pos_s + k ); // POS for DEL is 1-based pos before deletion
        fprintf(eva->vcf_out, ".\t"); // ID
        fprintf(eva->vcf_out, "%s\t", ref_seq); // REF
        fprintf(eva->vcf_out, "%s\t", path_seq); // ALT (should be the base before deletion if ALT is single base)
        fprintf(eva->vcf_out, ".\t"); // QUAL
        fprintf(eva->vcf_out, "PASS\t"); // FILTER
        strcpy(var_type, "DEL");
        snprintf(info_field, sizeof(info_field), "KMER_COV=%d;VARTYPE=%s", path_cov, var_type);
        fprintf(eva->vcf_out, "%s", info_field);
    }
    fprintf(eva->vcf_out, "\n");
    
    // fprintf(stdout, "##### ref_var.pos_s: %d \t", ref_var.pos_s);
    // fprintf(stdout, " ref_var.pos_t: %d \t", ref_var.pos_t);
    // fprintf(stdout, " new_var start_pos: %d \t", seq_var.pos_s);
    // fprintf(stdout, " new_var term_pos: %d \n\n\n\n\n", seq_var.pos_t);

    if (eva->fix_enabled) {
        if (eva->var_count >= eva->var_capacity) {
            eva->var_capacity = eva->var_capacity ? eva->var_capacity * 2 : 10;
            eva->variations = realloc(eva->variations, eva->var_capacity * sizeof(variation_t));
        }
        variation_t *var = &eva->variations[eva->var_count++];
        strncpy(var->chrom, eva->var_locs[var_loc_p].name, sizeof(var->chrom));
        var->chrom[sizeof(var->chrom)-1] = '\0';
        
        // Calculate genomic position based on k-mer start and k size
        var->pos = ref_var.pos_s + eva->k + 1; // 1-based end position of k-mer
        
        strncpy(var->ref, (const char*)ref_seq, sizeof(var->ref));
        var->ref[sizeof(var->ref)-1] = '\0';
        strncpy(var->alt, (const char*)path_seq, sizeof(var->alt));
        var->alt[sizeof(var->alt)-1] = '\0';
        
        if(slim_path_len > ref_seq_len) {
            strcpy(var->type, "INS");
        } else if(slim_path_len < ref_seq_len) {
            strcpy(var->type, "DEL");
        } else {
            strcpy(var->type, "SUB");
        }
    }
    return 0;
}

// Greedy path search function
int var_path_search_ref(evaluation_t *eva, int var_loc_p, int max_path_len, int assem_min_cov){ 
    debug_print("##### var_path_search_ref() \n");
    
    // Validate input parameters
    if (!eva || !eva->var_locs || var_loc_p < 0) {
        fprintf(stderr, "Error: Invalid parameters in var_path_search_ref\n");
        return 0;
    }
    
    kc_c4x_t *h = eva->h;
    int k = eva->k;
    int max_assem_cov = eva->max_assem_cov; // Get max coverage from eva struct
    uint64_t mask = (1ULL<<k*2) - 1;
    
    // Safely access var_locs
    var_location a_var = eva->var_locs[var_loc_p];
    uint64_t kmer_s = a_var.kmer_s;
    uint64_t kmer_t = a_var.kmer_t;
    int start_pos = a_var.pos_s;
    int term_pos = a_var.pos_t;
    
    // Validate positions
    if (start_pos < 0 || term_pos < 0 || start_pos >= term_pos) {
        fprintf(stderr, "Error: Invalid position values in var_path_search_ref: %d-%d\n", 
                start_pos, term_pos);
        return 0;
    }
    
    // Initialize the first node
    all_nodes[0].kmer = kmer_s;
    all_nodes[0].cov = kmer_cov(min_hash_key(kmer_s, k), mask, h);
    all_nodes[0].pre_node = NULL; // Ensure pre_node is initialized
    
    int a_cov = 1;
    if (a_cov < assem_min_cov) {
        a_cov = assem_min_cov;
    }
    
    term_nodes[0] = &all_nodes[0];
    
    int fresh_terms = 1, all_node_num = 1, p = 0;
    int good_term_node_num = 0;
    
    // Limit max_path_len to prevent buffer overflows
    if (max_path_len > 1000) {
        max_path_len = 1000;
        fprintf(stderr, "Warning: Limiting max_path_len to 1000\n");
    }
    
    path_node *p_node;
    path_node next_nodes[4];
    next_nodes[0].cov = 0; next_nodes[1].cov = 0; next_nodes[2].cov = 0; next_nodes[3].cov = 0;

    path_node *pre_term_nodes[Max_Path_Num+1]; 

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
                // Check both min and max coverage limits (max limit only if > 0)
                if(next_nodes[i].cov >=a_cov && (max_assem_cov <= 0 || next_nodes[i].cov <= max_assem_cov)){
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
int var_analysis_ref(evaluation_t *eva, int var_loc_p, int max_path_len, int min_cov, int assem_min_cov){ 
    debug_print("##### var_analysis_ref() \n");
    // debuging scripts
    // uint64_t mask = (1ULL<<eva->k*2) - 1;
    
    // Check for NULL pointers or invalid indices
    if (!eva || !eva->var_locs || !eva->kms || var_loc_p < 0) {
        fprintf(stderr, "Error: Invalid parameters in var_analysis_ref\n");
        return 0;
    }
    
    debug_print("\n\n\n\n\n#####\nVaration location: %d\n", var_loc_p);
    debug_print("\n%d\t", eva->var_locs[var_loc_p].pos_s);
    debug_print("%d\n", eva->var_locs[var_loc_p].pos_t);

    // Validate positions before proceeding
    if (eva->var_locs[var_loc_p].pos_s < 0 || eva->var_locs[var_loc_p].pos_t < 0) {
        fprintf(stderr, "Error: Invalid position values in var_locs[%d]\n", var_loc_p);
        return 0;
    }

    // greedy path search
    var_location ori_var;
    ori_var.kmer_s = eva->var_locs[var_loc_p].kmer_s;
    ori_var.kmer_t = eva->var_locs[var_loc_p].kmer_t;
    ori_var.pos_s = eva->var_locs[var_loc_p].pos_s;
    ori_var.pos_t = eva->var_locs[var_loc_p].pos_t;
    
    // Add bounds checking
    if (ori_var.pos_s >= ori_var.pos_t) {
        fprintf(stderr, "Error: Invalid position range in var_locs[%d]: %d-%d\n", 
                var_loc_p, ori_var.pos_s, ori_var.pos_t);
        return 0;
    }
    
    int good_term_node_num = var_path_search_ref(eva, var_loc_p, max_path_len, assem_min_cov);
    
    int ext_times = 0;
    while (good_term_node_num < 1 && ext_times < 5 && eva->var_locs[var_loc_p].pos_s > eva->k) {
        // Add safety check before extending
        if (eva->var_locs[var_loc_p].pos_s <= 0 || eva->var_locs[var_loc_p].pos_t <= 0) {
            fprintf(stderr, "Error: Invalid positions before extent_ var_loc\n");
            break;
        }
        
        // Use min_cov for extending variation locations
        extent_var_loc(eva, eva->kms, &eva->var_locs[var_loc_p], assem_min_cov, eva->k);
        
        // Verify positions after extension
        if (eva->var_locs[var_loc_p].pos_s < 0 || eva->var_locs[var_loc_p].pos_t < 0) {
            fprintf(stderr, "Error: Invalid positions after extent_ var_loc\n");
            break;
        }
        
        good_term_node_num = var_path_search_ref(eva, var_loc_p, max_path_len, assem_min_cov);
        ext_times = ext_times + 1;
    }
    
    if (good_term_node_num > 0) {
        // Validate good_term_node_num before using it
        if (good_term_node_num > Max_Path_Num) {
            fprintf(stderr, "Error: good_term_node_num (%d) exceeds Max_Path_Num (%d)\n", 
                    good_term_node_num, Max_Path_Num);
            return 0;
        }
        
        int best_termnode_index = best_term_node(eva, eva->var_locs[var_loc_p].pos_t - eva->var_locs[var_loc_p].pos_s, good_term_node_num);
        
        // Validate best_termnode_index
        if (best_termnode_index < 0 || best_termnode_index >= good_term_node_num) {
            fprintf(stderr, "Error: Invalid best_termnode_index: %d (max: %d)\n", 
                    best_termnode_index, good_term_node_num-1);
            return 0;
        }
        
        output_path(eva, var_loc_p, best_termnode_index);
        
        // if (good_term_node_num > 1) {
        //     fprintf(stdout, "\tMultiple-path");
        // }
        // fprintf(stdout, "\n");
    }else{
        // Output the location failed assembly to a separate file
        static FILE *low_quality_fp = NULL;
        static char low_quality_fn[256] = {0};
        
        // Open the low quality zones file if it's not already open
        if (low_quality_fp == NULL && eva->output_base != NULL) {
            snprintf(low_quality_fn, sizeof(low_quality_fn), "%s.low_quality_zones", eva->output_base);
            low_quality_fp = fopen(low_quality_fn, "w");
            if (!low_quality_fp) {
                fprintf(stderr, "Error: Could not open low quality zones file %s\n", low_quality_fn);
                return 0;
            }
            // Write header
            fprintf(low_quality_fp, "#CHROM\tSTART\tEND\tLENGTH\tITERATION\n");
        }
        
        if (low_quality_fp) {
            // Write location information to the file
            fprintf(low_quality_fp, "%s\t%d\t%d\t%d\t%d\n", 
                    eva->var_locs[var_loc_p].name,
                    eva->var_locs[var_loc_p].pos_s + eva->k + 1,  // Convert to 1-based genomic position
                    eva->var_locs[var_loc_p].pos_t + eva->k,      // End position
                    eva->var_locs[var_loc_p].pos_t - eva->var_locs[var_loc_p].pos_s,  // Length
                    eva->iteration);  // Current iteration
            
            // Flush the file to ensure data is written
            fflush(low_quality_fp);
        }
    }

    return good_term_node_num;
}

int best_term_node(evaluation_t *eva, int ref_node_num, int good_term_node_num){
    debug_print("\n ###best_term_node() \n");
    int i = 0, best_path=0, short_path_var_len = 10000, short_path_cov = 1;
    int p_path_node_num, var_len=0, p_path_cov;

    for(i=0; i < good_term_node_num; i++){
        p_path_node_num = path_node_num(good_term_nodes[i]);
        var_len = abs(p_path_node_num -ref_node_num) ;
        p_path_cov = nodes_path_cov(eva, good_term_nodes[i]);
        if(var_len < short_path_var_len || (p_path_node_num == short_path_var_len && p_path_cov > short_path_cov)){
            best_path = i;
            short_path_var_len = var_len;
            short_path_cov = p_path_cov;
        }
    }
    return best_path;
}

int extent_var_loc(evaluation_t *eva, uint64_t *kms, var_location *var_loc, int assem_min_cov, int dis){ 
    // Validate input parameters
    if (!eva || !kms || !var_loc) {
        fprintf(stderr, "Error: Invalid parameters\n");
        return 0;
    }
    
    int k = eva->k;
    uint64_t mask = (1ULL<<k*2) - 1;
    
    // Validate positions before extending
    if (var_loc->pos_s < 0 || var_loc->pos_t < 0) {
        fprintf(stderr, "Error: Invalid positions: %d-%d\n", 
                var_loc->pos_s, var_loc->pos_t);
        return 0;
    }
    
    // Add bounds checking to prevent going out of array bounds
    
    var_loc->pos_s = var_loc->pos_s - dis;
    if (var_loc->pos_s < 0) {
        var_loc->pos_s = 0;
    }
    
    var_loc->pos_t = var_loc->pos_t + dis;
    // We don't know the upper bound without knowing the size of kms array
    // So we'll rely on the cov check to prevent going too far
    
    int cov_s = kmer_cov(min_hash_key(kms[var_loc->pos_s], k), mask, eva->h);
    int cov_t = kmer_cov(min_hash_key(kms[var_loc->pos_t], k), mask, eva->h);
    
    int itt = 0;
    while (cov_s <= assem_min_cov && itt < dis && var_loc->pos_s > k) {
        var_loc->pos_s = var_loc->pos_s - 1;
        if (var_loc->pos_s < 0) {
            var_loc->pos_s = 0;
            break;
        }
        cov_s = kmer_cov(min_hash_key(kms[var_loc->pos_s], k), mask, eva->h);
        itt = itt + 1;
    }
    
    itt = 0;
    // We need to limit how far we can extend to avoid going out of bounds
    int max_extension = 1000; // Arbitrary limit
    while (cov_t <= assem_min_cov && itt < dis && itt < max_extension) {
        var_loc->pos_t = var_loc->pos_t + 1;
        cov_t = kmer_cov(min_hash_key(kms[var_loc->pos_t], k), mask, eva->h);
        itt = itt + 1;
    }
    
    var_loc->kmer_s = kms[var_loc->pos_s];
    var_loc->kmer_t = kms[var_loc->pos_t];
    
    // Check if we need to recursively extend
    int cov_rr = kmer_cov(min_hash_key(var_loc->kmer_s, k), mask, eva->hr);
    int cov_tt = kmer_cov(min_hash_key(var_loc->kmer_t, k), mask, eva->hr);
    
    // Limit recursion depth to prevent stack overflow
    static int recursion_depth = 0;
    if ((cov_rr > 1 || cov_tt > 1) && recursion_depth < 5) {
        recursion_depth++;
        extent_var_loc(eva, kms, var_loc, assem_min_cov, dis);
        recursion_depth--;
    }
    
    return 1;
} 


int optimize_var_loc(evaluation_t *eva, uint64_t *kms, var_location *var_loc, int assem_min_cov){ 
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
    while(cov_s > assem_min_cov && (cov_t/cov_s > 2.8 || (rep_km_num > 2 && shift_time <= k))){
        cov_s = kmer_cov(min_hash_key(kms[var_loc->pos_s - 1],k), mask, eva->h);
        if(cov_s >= assem_min_cov && cov_s > 0){
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

// Search for locations with variations
static evaluation_t *analysis_ref_seq(evaluation_t *eva, const char *fn, int max_len, int min_cov, int assem_min_cov)
{
    pldat_t pl;
    gzFile fp;
    int k = eva->k;

    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = eva->k;

    fprintf(stderr, "Analyzing the locations with variations ...... \n"); 
    uint64_t mask = (1ULL<<k*2) - 1;
    int ret;
    int i = 0;
    
    while ((ret = kseq_read(pl.ks)) >= 0) {// Reading one seq for each loop
        i = i + 1;
        fprintf(stderr, "Processing sequence %d: %s (length: %d)\n", i, pl.ks->name.s, pl.ks->seq.l);
        
        uint64_t *kms;  // all k-mers in array, all forward direction
        var_location *var_locs;
        
        int l = pl.ks->seq.l, km_num = 0; 
        if (l < k) {
            fprintf(stderr, "Skipping sequence shorter than k-mer size: %s\n", pl.ks->name.s);
            continue;
        }

        // Allocate memory for this sequence
        MALLOC(kms, l); 
        km_num = seq_kmers(kms, k, l, pl.ks->seq.s);
        
        MALLOC(var_locs, (l-k + 1));
        
        // Reset variation locations counter for this sequence
        int var_loc_num = 0;
        var_location one_var_loc;
        one_var_loc.kmer_s = 0;
        one_var_loc.pos_s = -1;
        
        // Find variation locations for this sequence
        int j = 0;
        while(j < km_num) {
            int cov = kmer_cov(min_hash_key(kms[j],k), mask, eva->h);
            
            if(cov < min_cov) { // k-mer NOT present 
                if(one_var_loc.pos_s < 0 && j > k) {
                    one_var_loc.kmer_s = kms[j-1];
                    one_var_loc.pos_s = j - 1;
                }
            } else { // k-mer present
                if(one_var_loc.pos_s > 0) {
                    int rep_km_num = kmer_cov(min_hash_key(kms[j],k), mask, eva->hr); 
                    int s_km_cov = kmer_cov(min_hash_key(one_var_loc.kmer_s,k), mask, eva->h);
                    
                    if(rep_km_num < 100 && s_km_cov/cov <= 100) {
                        var_locs[var_loc_num].kmer_s = one_var_loc.kmer_s;
                        var_locs[var_loc_num].pos_s = one_var_loc.pos_s;
                        var_locs[var_loc_num].kmer_t = kms[j];
                        var_locs[var_loc_num].pos_t = j;
                        
                        strncpy(var_locs[var_loc_num].name, pl.ks->name.s, sizeof(var_locs[var_loc_num].name)-1);
                        var_locs[var_loc_num].name[sizeof(var_locs[var_loc_num].name)-1] = '\0';
                        
                        if(var_locs[var_loc_num].pos_t - var_locs[var_loc_num].pos_s <= k) {
                            var_locs[var_loc_num].pos_t = var_locs[var_loc_num].pos_s + k + 1;
                            var_locs[var_loc_num].kmer_t = kms[var_locs[var_loc_num].pos_t];
                        }
                        
                        one_var_loc.kmer_s = 0;
                        one_var_loc.pos_s = -1;
                        var_loc_num++;
                    }
                }
            }
            j++;
        }
        
        // Store the data for this sequence in the evaluation structure
        eva->var_locs = var_locs;
        eva->kms = kms;
        
        fprintf(stderr, "Found %d variation locations in sequence: %s\n", var_loc_num, pl.ks->name.s);
        


        // Process variations for this sequence
        fprintf(stderr, "Analyzing the variation details for sequence: %s\n", pl.ks->name.s);
        int var_loc_p = 0;
        while(var_loc_p < var_loc_num) { 
            // Pass both min_cov and assem_min_cov
            var_analysis_ref(eva, var_loc_p, max_len, min_cov, assem_min_cov); 
            var_loc_p++;
        }
        
        // Free memory for this sequence before moving to the next one
        free(kms);
        free(var_locs);
        eva->kms = NULL;
        eva->var_locs = NULL;
    }

    kseq_destroy(pl.ks);
    gzclose(fp);
    
    return eva;
}

// Complete the apply_variations function to include INS and DEL handling
void apply_variations(evaluation_t *eva, const char *ref_file, const char *output_base) {
    char *output_fn = strdup(output_base);
    
    gzFile fp = gzopen(ref_file, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open reference file %s\n", ref_file);
        free(output_fn);
        return;
    }

    kseq_t *seq = kseq_init(fp);
    FILE *out_fp = fopen(output_fn, "w");
    if (!out_fp) {
        fprintf(stderr, "Failed to open output file %s\n", output_fn);
        free(output_fn);
        return;
    }

    // Sort variations by chromosome (sequence name) and position
    qsort(eva->variations, eva->var_count, sizeof(variation_t), compare_variations);
  
    // No longer tracking curr_var_idx; removed unused variable
    
    while (kseq_read(seq) >= 0) {
        char *seqname = seq->name.s;
        char *sequence = strdup(seq->seq.s);
        int len = seq->seq.l;
        
        fprintf(stderr, "Processing sequence: %s (length: %d)\n", seqname, len);
        
        // Find variations that apply to this sequence
        int seq_var_count = 0;
        variation_t *seq_variations = NULL;
        
        for (int i = 0; i < eva->var_count; i++) {
            if (strcmp(eva->variations[i].chrom, seqname) == 0) {
                seq_var_count++;
            }
        }
        
        if (seq_var_count > 0) {
            // Allocate memory for sequence-specific variations
            MALLOC(seq_variations, seq_var_count);
            int j = 0;
            
            for (int i = 0; i < eva->var_count; i++) {
                if (strcmp(eva->variations[i].chrom, seqname) == 0) {
                    seq_variations[j++] = eva->variations[i];
                }
            }
            
            // Sort sequence-specific variations by position (descending)
            qsort(seq_variations, seq_var_count, sizeof(variation_t), compare_variations);
            
            // Apply variations to this sequence
            for (int i = 0; i < seq_var_count; i++) {
                variation_t *var = &seq_variations[i];
                int pos = var->pos - 1; // Convert to 0-based
                
                if (pos < 0 || pos >= len) {
                    fprintf(stderr, "  Skip: Position %d out of bounds for sequence %s (length: %d)\n", 
                            var->pos, seqname, len);
                    continue;
                }

                // Debug output before applying variation
                // fprintf(stderr, "Applying variation at position %d:\n", var->pos);
                // fprintf(stderr, "  Type: %s\n", var->type);
                // fprintf(stderr, "  REF: '%s'\n", var->ref);
                // fprintf(stderr, "  ALT: '%s'\n", var->alt);
                // fprintf(stderr, "  Current sequence at position: '%.10s'\n", &sequence[pos]);

                if (strcmp(var->type, "SUB") == 0) {
                    int ref_len = strlen(var->ref);
                    int alt_len = strlen(var->alt);
                    if (pos + ref_len > len) {
                        fprintf(stderr, "  Skip: Reference length exceeds sequence bounds\n");
                        continue;
                    }
                    if (strncmp(&sequence[pos], var->ref, ref_len) != 0) {
                        fprintf(stderr, "  REF mismatch at %s:%d\n", seqname, var->pos);
                        fprintf(stderr, "    Expected: '%s'\n", var->ref);
                        fprintf(stderr, "    Found: '%.10s'\n", &sequence[pos]);
                        continue;
                    }
                    // Replace REF with ALT
                    memcpy(&sequence[pos], var->alt, alt_len);
                    // fprintf(stderr, "  After substitution: '%.10s'\n", &sequence[pos]);
                } else if (strcmp(var->type, "INS") == 0) {
                    // Insert ALT after POS (VCF format)
                    int ref_len = strlen(var->ref);
                    int alt_len = strlen(var->alt);
                    int ins_len = alt_len - ref_len;  // True insertion length
                    
                    // Check for REF mismatch if REF has content
                    if (ref_len > 0) {
                        if (pos + ref_len > len) {
                            fprintf(stderr, "  Skip: Reference length exceeds sequence bounds\n");
                            continue;
                        }
                        if (strncmp(&sequence[pos], var->ref, ref_len) != 0) {
                            fprintf(stderr, "  REF mismatch at %s:%d\n", seqname, var->pos);
                            fprintf(stderr, "    Expected: '%s'\n", var->ref);
                            fprintf(stderr, "    Found: '%.10s'\n", &sequence[pos]);
                            continue;
                        }
                    }
                    
                    // For insertions, we need to allocate more space
                    char *new_seq = malloc(len + ins_len + 1);
                    if (!new_seq) {
                        fprintf(stderr, "  Error: Memory allocation failed for insertion\n");
                        continue;
                    }
                    
                    // Copy sequence up to the insertion point
                    memcpy(new_seq, sequence, pos + ref_len);
                    // Insert ALT (which includes the REF part)
                    memcpy(new_seq + pos, var->alt, alt_len);
                    // Copy the remainder of the sequence
                    memcpy(new_seq + pos + alt_len, sequence + pos + ref_len, len - pos - ref_len);
                    // Null terminate
                    new_seq[len + ins_len] = '\0';
                    
                    // Clean up and update
                    free(sequence);
                    sequence = new_seq;
                    len += ins_len;
                    
                    // fprintf(stderr, "  After insertion: '%.10s'\n", &sequence[pos]);
                } else if (strcmp(var->type, "DEL") == 0) {
                    int ref_len = strlen(var->ref);
                    int alt_len = strlen(var->alt);
                    int del_len = ref_len - alt_len;
                    
                    if (pos + ref_len > len || del_len <= 0) {
                        fprintf(stderr, "  Skip: Invalid deletion parameters\n");
                        continue;
                    }
                    
                    if (strncmp(&sequence[pos], var->ref, ref_len) != 0) {
                        fprintf(stderr, "  REF mismatch at %s:%d\n", seqname, var->pos);
                        fprintf(stderr, "    Expected: '%s'\n", var->ref);
                        fprintf(stderr, "    Found: '%.10s'\n", &sequence[pos]);
                        continue;
                    }
                    
                    // Replace with ALT and delete remaining bases
                    memcpy(&sequence[pos], var->alt, alt_len);
                    memmove(&sequence[pos + alt_len], &sequence[pos + ref_len], len - pos - ref_len + 1); // +1 for null terminator
                    len -= del_len;
                    
                    // fprintf(stderr, "  After deletion: '%.10s'\n", &sequence[pos]);
                }
            }
            
            // Free sequence-specific variations
            free(seq_variations);
        }

        // Write modified sequence to output
        fprintf(out_fp, ">%s\n%s\n", seq->name.s, sequence);
        free(sequence);
    }

    kseq_destroy(seq);
    gzclose(fp);
    fclose(out_fp);
    fprintf(stderr, "Successfully wrote corrected genome to %s\n", output_fn);
    free(output_fn);
}





// Fix memory management in main iteration loop
int main(int argc, char *argv[])
{
    int c, k = 31, p = KC_BITS, block_size = 10000000, n_thread = 3, min_cov = 2, assem_min_cov = 5, insert_size = 1000;
    float error_rate = 0.025f;
    int num_iters = 3;  // Default number of iterations
    int max_assem_cov = 0; // Default: no max coverage limit
    char *output_base = NULL;

    ketopt_t o = KETOPT_INIT;
    int fix_enabled = 0;
    static ko_longopt_t long_options[] = {
        { "fix", ko_optional_argument, 128 },
        { 0, 0, 0 }
    };
    
    while ((c = ketopt(&o, argc, argv, 1, "k:t:c:a:l:e:o:n:m:", long_options)) >= 0) { // Added 'm:'
        if (c == 'k') k = atoi(o.arg);
        else if (c == 't') n_thread = atoi(o.arg);
        else if (c == 'c') min_cov = atoi(o.arg);
        else if (c == 'a') assem_min_cov = atoi(o.arg);  // New option for assembly min coverage
        else if (c == 'm') max_assem_cov = atoi(o.arg); // Parse max assembly coverage
        else if (c == 'l') insert_size = atoi(o.arg); 
        else if (c == 'e') error_rate = atof(o.arg);
        else if (c == 'n') num_iters = atoi(o.arg);  // New option for iterations
        else if (c == 128) {
            fix_enabled = 1;
        }
        else if (c == 'o') output_base = o.arg;
    }

    if (argc - o.ind < 2) {
        usage(k, n_thread, min_cov, assem_min_cov, insert_size, error_rate, max_assem_cov); // Pass max_assem_cov
        return 1;
    }

    if (!output_base) {
        fprintf(stderr, "Error: Output base name (-o) is required\n");
        return 1;
    }

    fprintf(stderr, "k-mer size: %d\n", k);
    
    // Calculate the number of files
    int f_num = argc - o.ind;
    
    // Check if reference file exists
    if (access(argv[o.ind], F_OK) != 0) {
        fprintf(stderr, "Error: Reference file '%s' does not exist or is not accessible\n", argv[o.ind]);
        return 1;
    }

    // Check if all NGS files exist before starting
    for (int i = 1; i < f_num; i++) {
        if (access(argv[o.ind + i], F_OK) != 0) {
            fprintf(stderr, "Error: NGS file '%s' does not exist or is not accessible\n", argv[o.ind + i]);
            return 1;
        }
    }
    
    fprintf(stderr, "Counting k-mers of NGS file 1  ......\n");
    
    kc_c4x_t *h, *hr;
    h = count_file(argv[o.ind + 1], k, p, block_size, n_thread);
    
    int c_f_n = 2;
    while(c_f_n < f_num ) {
        fprintf(stderr, "Counting NGS file %d ......\n", c_f_n );
        h = count_file2(argv[o.ind + c_f_n ], h, k, p, block_size, n_thread);  
        c_f_n = c_f_n + 1;
    }

    // Update the print_hist call to include the output_base parameter
    print_hist(h, n_thread, output_base);

    //fprintf(stderr, "Analyzing k-mers of reference sequence ......\n");
    //hr = count_file(argv[o.ind + 1], k, p, block_size, n_thread);

    fprintf(stderr, "Analyzing k-mers of reference sequence ......\n");
    hr = count_file(argv[o.ind], k, p, block_size, n_thread);

    evaluation_t eva;
    eva.h = h;
    eva.hr = hr;
    eva.k = k;
    eva.p = p;
    eva.error_rate = error_rate;
    eva.fix_enabled = fix_enabled;
    eva.max_assem_cov = max_assem_cov; // Set max coverage in eva struct
    eva.var_count = 0;
    eva.var_capacity = 0;
    eva.variations = NULL;
    eva.output_base = output_base;  // Set the output_base in the eva struct

    char *current_ref = argv[o.ind];
    for (int iter = 1; iter <= num_iters; iter++) {
        // Open VCF output for this iteration
        char vcf_fn[256];
        snprintf(vcf_fn, sizeof(vcf_fn), "%s.iter%d.vcf", output_base, iter);
        FILE *vcf_fp = fopen(vcf_fn, "w");
        if (!vcf_fp) {
            fprintf(stderr, "Error: Could not open VCF file %s\n", vcf_fn);
            continue;
        }
        
        fprintf(vcf_fp, "##fileformat=VCFv4.5\n");
        fprintf(vcf_fp, "##ProGenFixerVersion=v1.0\n");
        fprintf(vcf_fp, "##ProGenFixerCommand=");
        for (int i = 0; i < argc; i++) {
            fprintf(vcf_fp, "%s%s", argv[i], (i == argc - 1) ? "" : " ");
        }
        fprintf(vcf_fp, "\n");

        fprintf(vcf_fp, "##INFO=<ID=KMER_COV,Number=1,Type=Integer,Description=\"K-mer coverage of the variant path\">\n");
        fprintf(vcf_fp, "##INFO=<ID=VARTYPE,Number=1,Type=String,Description=\"Variant type (INS, DEL, SUB)\">\n");
        fprintf(vcf_fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

        // Count current reference k-mers
        fprintf(stderr, "Counting k-mers in reference for iteration %d...\n", iter);
        kc_c4x_t *hr = count_file(current_ref, k, p, block_size, n_thread);
        if (!hr) {
            fprintf(stderr, "Error: Failed to count k-mers in reference %s\n", current_ref);
            fclose(vcf_fp);
            continue;
        }

        // Initialize evaluation struct
        eva.hr = hr;
        eva.vcf_out = vcf_fp;
        eva.iteration = iter;
        eva.var_count = 0;
        eva.var_capacity = 0;
        
        // Free previous variations if any
        if (eva.variations != NULL) {
            free(eva.variations);
        }
        eva.variations = NULL;
        
        // Reset pointers to prevent dangling pointers
        eva.var_locs = NULL;
        eva.kms = NULL;

        // Run analysis
        fprintf(stderr, "Running analysis for iteration %d...\n", iter);
        analysis_ref_seq(&eva, current_ref, insert_size, min_cov, assem_min_cov);
        fclose(vcf_fp);

        // Apply corrections if we found any variations
        if (eva.var_count > 0) {
            char new_ref[256];
            snprintf(new_ref, sizeof(new_ref), "%s.iter%d.fasta", output_base, iter);
            fprintf(stderr, "Applying %d variations to create %s\n", eva.var_count, new_ref);
            apply_variations(&eva, current_ref, new_ref);
            
            // Update current_ref for next iteration
            if (iter > 1) {
                free(current_ref);
            }
            current_ref = strdup(new_ref);
            fprintf(stderr, "Next reference: %s\n", current_ref);
        } else {
            fprintf(stderr, "No variations found in iteration %d, using same reference for next iteration\n", iter);
        }
        
        // Cleanup from this iteration
        if (eva.variations != NULL) {
            free(eva.variations);
            eva.variations = NULL;
        }
        
        // Free the reference hash for this iteration
        for (int i = 0; i < (1<<p); i++) {
            kc_c4_destroy(hr->h[i]);
        }
        free(hr->h);
        free(hr);
    }

    return 0;
}

void usage(int k, int n_thread, int min_cov, int assem_min_cov, int insert_size, float error_rate, int max_assem_cov) { // Added max_assem_cov parameter
    fprintf(stderr, "\n");
    fprintf(stderr, "ProGenFixer: an ultra-fast and accurate tool for correcting prokaryotic genome sequences using a mapping-free algorithm\n");
    fprintf(stderr, "Version v1.0  Author: Lifu Song songlf@tib.cas.cn\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: ProGenFixer [options] Reference NGS_files -o output_prefix\n");
    fprintf(stderr, "       [Supporting formats: fq, fa, fq.gz, fa.gz]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
    fprintf(stderr, "  -c INT     minimal k-mer coverage for identifying variation regions [%d]\n", min_cov);
    fprintf(stderr, "  -a INT     minimal k-mer coverage for assemble-based variant calling [%d]\n", assem_min_cov);
    fprintf(stderr, "  -m INT     maximal k-mer coverage for assemble-based variant calling [%d, 0=no limit]\n", max_assem_cov); // Added -m description
    fprintf(stderr, "  -l INT     maximal assembly length [%d]\n", insert_size);
    fprintf(stderr, "  -t INT     number of threads [%d]\n", n_thread);
    // fprintf(stderr, "  -e FLOAT   sequencing error rate for p-value calculation [%g]\n", error_rate);
    fprintf(stderr, "  -n INT     number of correction iterations [3]\n");
    fprintf(stderr, "  -o STR     base name for output files [required]\n");
    fprintf(stderr, "  --fix      enable reference correction \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  ProGenFixer ref_genome.fa reads.fq -o results --fix\n");
    fprintf(stderr, "  ProGenFixer ref_genome.fa reads1.fq reads2.fq -o results -a 5 --fix\n");
    fprintf(stderr, "  ProGenFixer ref_genome.fa reads1.fq reads2.fq -o results -c 4 -a 6 -m 100 -n 3 --fix\n"); // Added example with -m
    fprintf(stderr, "\n");
}



