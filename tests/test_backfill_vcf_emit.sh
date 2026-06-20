#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/progenfixer-backfill.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT

cat >"$tmpdir/backfill_emit_harness.c" <<C
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define main progenfixer_main
#include "$repo_root/ProGenFixer.c"
#undef main

static int encode_one(const char *seq, int k, uint64_t *out)
{
    uint64_t kms[8];
    int n = seq_kmers(kms, k, k, seq);
    if (n != 1) return -1;
    *out = kms[0];
    return 0;
}

int main(int argc, char **argv)
{
    if (argc != 2) return 2;

    char ref_path[4096];
    char vcf_path[4096];
    char branch_ref_path[4096];
    char branch_reads_path[4096];
    snprintf(ref_path, sizeof(ref_path), "%s/ref.fa", argv[1]);
    snprintf(vcf_path, sizeof(vcf_path), "%s/backfill.vcf", argv[1]);
    snprintf(branch_ref_path, sizeof(branch_ref_path), "%s/branch_ref.fa", argv[1]);
    snprintf(branch_reads_path, sizeof(branch_reads_path), "%s/branch_reads.fa", argv[1]);

    FILE *ref = fopen(ref_path, "w");
    if (ref == NULL) return 3;
    fprintf(ref, ">ref\\nAAACCCGTA\\n");
    fclose(ref);

    ref_pos_ctx_t ctx;
    if (ref_pos_ctx_build(ref_path, 3, &ctx) != 0) return 4;

    FILE *vcf = fopen(vcf_path, "w");
    if (vcf == NULL) return 5;

    evaluation_t eva;
    memset(&eva, 0, sizeof(eva));
    eva.k = 3;
    eva.vcf_out = vcf;

    uint64_t left_fwd = 0, right_fwd = 0;
    if (encode_one("AAA", 3, &left_fwd) != 0 ||
        encode_one("CCC", 3, &right_fwd) != 0) {
        return 6;
    }

    const unsigned char assembled[] = "AAAGGCCC";
    int emitted = emit_backfill_vcf(&eva, &ctx,
                                    min_hash_key(left_fwd, 3), left_fwd,
                                    min_hash_key(right_fwd, 3), right_fwd,
                                    assembled, (int)strlen((const char *)assembled),
                                    7);
    fclose(vcf);
    ref_pos_ctx_free(&ctx);

    FILE *branch_ref = fopen(branch_ref_path, "w");
    if (branch_ref == NULL) return 8;
    fprintf(branch_ref, ">ref\\nAAACCC\\n");
    fclose(branch_ref);

    FILE *branch_reads = fopen(branch_reads_path, "w");
    if (branch_reads == NULL) return 9;
    for (int i = 0; i < 20; i++) fprintf(branch_reads, ">dead_%d\\nAAAGTAC\\n", i);
    for (int i = 0; i < 5; i++) fprintf(branch_reads, ">real_%d\\nAAAGGCCC\\n", i);
    fclose(branch_reads);

    evaluation_t branch_eva;
    memset(&branch_eva, 0, sizeof(branch_eva));
    branch_eva.k = 3;
    branch_eva.h = count_file(branch_reads_path, 3, KC_BITS, 100000, 1);
    branch_eva.hr = count_file(branch_ref_path, 3, KC_BITS, 100000, 1);
    branch_eva.used_kmers = u64set_init();

    uint64_t start_fwd = 0, expected_anchor = 0;
    uint64_t expected_path[2], found_path[16], found_anchor = 0;
    int found_len = 0;
    if (encode_one("AAG", 3, &start_fwd) != 0 ||
        encode_one("CCC", 3, &expected_anchor) != 0 ||
        encode_one("AGG", 3, &expected_path[0]) != 0 ||
        encode_one("GGG", 3, &expected_path[1]) != 0) {
        return 10;
    }

    int found = backfill_extend(&branch_eva, start_fwd, +1, 5, 2,
                                found_path, &found_len, &found_anchor);
    if (found != 1 || found_len != 2 ||
        found_anchor != min_hash_key(expected_anchor, 3)) {
        return 11;
    }
    for (int i = 0; i < 2; i++) {
        if (found_path[i] != expected_path[i]) return 12;
    }

    free_variations(&eva);

    return emitted == 1 ? 0 : 7;
}
C

cc -g -Wall -O0 -I"$repo_root" \
    -o "$tmpdir/backfill_emit_harness" \
    "$tmpdir/backfill_emit_harness.c" "$repo_root/kthread.c" \
    -lz -lm -lpthread

"$tmpdir/backfill_emit_harness" "$tmpdir"

expected=$'ref\t3\t.\tA\tAGG\t.\tPASS\tKMER_COV=7;VARTYPE=INS;SOURCE=BACKFILL'
if ! grep -Fxq "$expected" "$tmpdir/backfill.vcf"; then
    cat "$tmpdir/backfill.vcf" >&2
    exit 1
fi

cat >"$tmpdir/ref_cli.fa" <<'EOF'
>ref
AAACCCGTA
EOF

{
    for i in 1 2 3 4 5; do printf ">ref_%s\nAAACCCGTA\n" "$i"; done
    for i in 1 2 3 4 5; do printf ">alt_%s\nAAAGGCCCGTA\n" "$i"; done
} >"$tmpdir/reads_cli.fa"

make -C "$repo_root" >/dev/null

"$repo_root/ProGenFixer" \
    -k 3 -c 3 -a 5 -n 1 \
    -o "$tmpdir/out_cli" \
    "$tmpdir/ref_cli.fa" "$tmpdir/reads_cli.fa" \
    >"$tmpdir/cli.stdout" \
    2>"$tmpdir/cli.stderr"

if ! grep -q "Found 0 variation locations" "$tmpdir/cli.stderr"; then
    cat "$tmpdir/cli.stderr" >&2
    exit 1
fi

if ! grep -q "Back-fill (iter 1): unique-NGS=3, novel=1, contamination=0, vcf-emitted=1" "$tmpdir/cli.stderr"; then
    cat "$tmpdir/cli.stderr" >&2
    exit 1
fi

expected_cli=$'ref\t3\t.\tA\tAGG\t.\tPASS\tKMER_COV=5;VARTYPE=INS;SOURCE=BACKFILL'
if ! grep -Fxq "$expected_cli" "$tmpdir/out_cli.iter1.vcf"; then
    cat "$tmpdir/out_cli.iter1.vcf" >&2
    exit 1
fi

echo "backfill VCF emit regression passed"
