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
    snprintf(ref_path, sizeof(ref_path), "%s/ref.fa", argv[1]);
    snprintf(vcf_path, sizeof(vcf_path), "%s/backfill.vcf", argv[1]);

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
