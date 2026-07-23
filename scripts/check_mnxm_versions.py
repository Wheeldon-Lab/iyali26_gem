"""Run this in the project root to verify the MNXM-version mismatch hypothesis.

Hypothesis: chem_prop.tsv and reac_prop.tsv reference different MNXM
namespaces (MetaNetX 4.4 7-digit IDs vs older 5-digit IDs), so the
fingerprint_index built from reac_prop.tsv equations cannot be matched
against the MNXM IDs annotated on metabolites (which came from chem_prop.tsv
or chem_xref.tsv).
"""
from pathlib import Path
import re
from collections import Counter

MNX_DIR = Path("data/metanetx")  # adjust if different

def mnxm_lengths_in_file(path, sample_n=20000):
    """Return Counter of MNXM ID digit-lengths found in file."""
    pat = re.compile(r"MNXM(\d+)")
    counter = Counter()
    examples_by_len = {}
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= sample_n:
                break
            for m in pat.finditer(line):
                digits = m.group(1)
                counter[len(digits)] += 1
                if len(digits) not in examples_by_len:
                    examples_by_len[len(digits)] = m.group(0)
    return counter, examples_by_len

print("=== MNXM ID digit-length distribution ===\n")

for fname in ("chem_prop.tsv", "chem_xref.tsv", "reac_prop.tsv", "reac_xref.tsv"):
    path = MNX_DIR / fname
    if not path.exists():
        print(f"{fname}: NOT FOUND at {path}")
        continue
    cnt, ex = mnxm_lengths_in_file(path)
    print(f"{fname}:")
    for length in sorted(cnt):
        print(f"  {length}-digit MNXM: {cnt[length]:>10,d}  example: {ex[length]}")
    print()

# Now check: pick a metabolite from the model that DID get an MNXM assigned
# (e.g. m184 → MNXM1093693), and see if that exact MNXM ID appears in
# reac_prop.tsv.
print("=== Spot-check: do model MNXM IDs appear in reac_prop.tsv? ===\n")
test_ids = ["MNXM1093693", "MNXM735438", "MNXM732620", "MNXM1108357",
            # Old-namespace IDs that MIGHT be the "real" ones in reac_prop
            "MNXM48", "MNXM4", "MNXM2", "MNXM22"]

reac_prop = MNX_DIR / "reac_prop.tsv"
if reac_prop.exists():
    text = reac_prop.read_text()
    for tid in test_ids:
        # Match as whole token (followed by @ in equation, or end of word)
        hits = len(re.findall(rf"\b{tid}\b", text))
        print(f"  {tid}: {hits} occurrence(s) in reac_prop.tsv")
else:
    print(f"  {reac_prop} not found — adjust MNX_DIR path at top of script")
