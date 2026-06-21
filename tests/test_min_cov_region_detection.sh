#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/progenfixer-mincov.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT

python3 - "$tmpdir" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
k = 5
random.seed(7)
comp = str.maketrans("ACGT", "TGCA")

def canonical(seq):
    rc = seq.translate(comp)[::-1]
    return min(seq, rc)

seq = "".join(random.choice("ACGT") for _ in range(k))
seen = {canonical(seq)}
while len(seq) < 80:
    choices = []
    for base in "ACGT":
        mer = seq[-(k - 1):] + base
        can = canonical(mer)
        if can not in seen:
            choices.append((base, can))
    if not choices:
        raise RuntimeError("failed to build unique reference sequence")
    base, can = random.choice(choices)
    seq += base
    seen.add(can)

kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
coverage = [5] * len(kmers)

# Two true low-coverage blocks separated by weak coverage. The weak middle
# block must not become an artificial region boundary when -c is low.
for i in range(20, 26):
    coverage[i] = 1
for i in range(26, 32):
    coverage[i] = 3
for i in range(32, 38):
    coverage[i] = 1

# A separate weak block should become a candidate only when -c is high enough.
for i in range(45, 51):
    coverage[i] = 3

(out / "ref.fa").write_text(">ref\n" + seq + "\n")
with (out / "reads.fa").open("w") as fh:
    read_id = 0
    for idx, mer in enumerate(kmers):
        for _ in range(coverage[idx]):
            fh.write(f">r{read_id}_k{idx}\n{mer}\n")
            read_id += 1
PY

make -C "$repo_root" >/dev/null

extract_locations() {
    awk '/Found [0-9]+ variation locations/ { print $2; exit }' "$1"
}

for min_cov in 2 4; do
    "$repo_root/ProGenFixer" \
        -k 5 -a 5 -c "$min_cov" -n 1 --no-backfill \
        -o "$tmpdir/out_c${min_cov}" \
        "$tmpdir/ref.fa" "$tmpdir/reads.fa" \
        >"$tmpdir/c${min_cov}.stdout" \
        2>"$tmpdir/c${min_cov}.stderr"
done

count_c2="$(extract_locations "$tmpdir/c2.stderr")"
count_c4="$(extract_locations "$tmpdir/c4.stderr")"

if [[ "$count_c2" != "1" ]]; then
    echo "expected -c 2 to find 1 variation location, got $count_c2" >&2
    exit 1
fi

if [[ "$count_c4" != "2" ]]; then
    echo "expected -c 4 to find 2 variation locations, got $count_c4" >&2
    exit 1
fi

echo "min_cov region detection regression passed"
