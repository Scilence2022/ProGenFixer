#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/progenfixer-apply.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT

cat >"$tmpdir/apply_variations_harness.c" <<C
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define main progenfixer_main
#include "$repo_root/ProGenFixer.c"
#undef main

static void write_repeat(FILE *fp, char base, int n)
{
    for (int i = 0; i < n; i++) fputc(base, fp);
}

int main(int argc, char **argv)
{
    if (argc != 2) return 2;

    char ref_path[4096];
    char out_path[4096];
    snprintf(ref_path, sizeof(ref_path), "%s/ref.fa", argv[1]);
    snprintf(out_path, sizeof(out_path), "%s/fixed.fa", argv[1]);

    char long_name[180];
    memset(long_name, 'L', sizeof(long_name) - 1);
    long_name[sizeof(long_name) - 1] = '\\0';

    FILE *fp = fopen(ref_path, "w");
    if (!fp) return 3;
    fprintf(fp, ">%s\\nAACCGGTT\\n", long_name);
    fprintf(fp, ">big\\n");
    write_repeat(fp, 'A', 1100);
    fprintf(fp, "TT\\n");
    fclose(fp);

    evaluation_t eva;
    memset(&eva, 0, sizeof(eva));
    eva.fix_enabled = 1;

    var_location ref_var;
    memset(&ref_var, 0, sizeof(ref_var));
    ref_var.pos_s = 257867;
    if (varcall_vcf_pos(&ref_var, 31, "DEL") != 257899) return 5;
    if (varcall_vcf_pos(&ref_var, 31, "SUB") != 257899) return 6;

    record_variation(&eva, long_name, 3, "C", 1, "T", 1, "SUB");
    record_variation(&eva, long_name, 3, "C", 1, "T", 1, "SUB");

    char *big_ref = malloc(1101);
    char *big_alt = malloc(1103);
    if (!big_ref || !big_alt) return 4;
    memset(big_ref, 'A', 1100);
    big_ref[1100] = '\\0';
    memset(big_alt, 'A', 1100);
    big_alt[1100] = 'C';
    big_alt[1101] = 'C';
    big_alt[1102] = '\\0';
    record_variation(&eva, "big", 1, big_ref, 1100, big_alt, 1102, "INS");
    free(big_ref);
    free(big_alt);

    apply_variations(&eva, ref_path, out_path);
    free_variations(&eva);
    return 0;
}
C

cc -g -Wall -O0 -I"$repo_root" \
    -o "$tmpdir/apply_variations_harness" \
    "$tmpdir/apply_variations_harness.c" "$repo_root/kthread.c" \
    -lz -lm -lpthread

"$tmpdir/apply_variations_harness" "$tmpdir" >/dev/null 2>"$tmpdir/harness.stderr"

if grep -q "REF mismatch" "$tmpdir/harness.stderr"; then
    cat "$tmpdir/harness.stderr" >&2
    exit 1
fi

python3 - "$tmpdir/fixed.fa" <<'PY'
import sys
from pathlib import Path

records = {}
name = None
parts = []
for line in Path(sys.argv[1]).read_text().splitlines():
    if line.startswith(">"):
        if name is not None:
            records[name] = "".join(parts)
        name = line[1:]
        parts = []
    else:
        parts.append(line.strip())
if name is not None:
    records[name] = "".join(parts)

long_name = "L" * 179
assert records[long_name] == "AATCGGTT", records.get(long_name)
assert records["big"] == ("A" * 1100) + "CCTT", len(records["big"])
PY

echo "apply_variations regression passed"
