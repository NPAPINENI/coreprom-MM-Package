"""
Unified Pipeline
================
Inputs  : <region>  <infile.tsv>  <outdir>
           - infile must have columns: Sequence, Label

Steps   :
  1. Uppercase + validate sequences
  2. Generate 6-mers (sliding window = 1)
  3. Stratified 70 / 15 / 15 split  →  sequence_train/test/dev.tsv
  4. Compute sequence composition features  (A/C/G/T fractions, CpG1/2/3, etc.)
  5. Compute DNAshapeR MAD features for 14 shape parameters  (via rpy2)
  6. Write properties_train/test/dev.tsv  (row-aligned with sequence files)

Usage   :
  python pipeline.py enhancer input.tsv ./output/

Dependencies (Python):
  pip install pandas scikit-learn numpy rpy2

Dependencies (R — auto-installed if missing):
  BiocManager, Biostrings, DNAshapeR
"""

import sys
import os
import tempfile
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

# ── rpy2 imports ──────────────────────────────────────────────────────────────
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects.vectors import StrVector
    pandas2ri.activate()
    RPY2_AVAILABLE = True
except ImportError:
    RPY2_AVAILABLE = False
    print("⚠️  rpy2 not found — DNA shape features will be skipped.")
    print("    Install with: pip install rpy2")

# ─────────────────────────────────────────────────────────────────────────────
K           = 6
TRAIN_RATIO = 0.70
TEST_RATIO  = 0.15
DEV_RATIO   = 0.15
RANDOM_SEED = 42

SHAPE_FEATURES = [
    "MGW", "Roll", "ProT", "HelT", "EP",
    "Rise", "Shift", "Slide", "Tilt", "Buckle",
    "Opening", "Shear", "Stagger", "Stretch"
]
# ─────────────────────────────────────────────────────────────────────────────


# ══════════════════════════════════════════════════════════════════════════════
# 1 – K-MER GENERATION
# ══════════════════════════════════════════════════════════════════════════════

def generate_kmers(sequence: str, k: int = 6) -> str:
    """Return space-joined k-mers with sliding window of 1."""
    return " ".join(sequence[i:i + k] for i in range(len(sequence) - k + 1))


# ══════════════════════════════════════════════════════════════════════════════
# 2 – SEQUENCE COMPOSITION FEATURES
# ══════════════════════════════════════════════════════════════════════════════

def compute_composition_features(sequence: str, label: int,
                                 kmers: str, region: str) -> dict:
    """
    Nucleotide composition features computed on the FULL sequence.
    Returns a flat dict (one row).
    """
    seq     = sequence.upper()
    n       = len(seq)

    A = seq.count('A'); C = seq.count('C')
    G = seq.count('G'); T = seq.count('T')

    GC = seq.count("GC"); CG = seq.count("CG")

    ACG = seq.count("ACG"); AGC = seq.count("AGC"); CAG = seq.count("CAG")
    CCG = seq.count("CCG"); CGA = seq.count("CGA"); CGC = seq.count("CGC")
    CGG = seq.count("CGG"); CGT = seq.count("CGT"); CTG = seq.count("CTG")
    GAC = seq.count("GAC"); GCA = seq.count("GCA"); GCC = seq.count("GCC")
    GCG = seq.count("GCG"); GCT = seq.count("GCT"); GGC = seq.count("GGC")
    GTC = seq.count("GTC"); TCG = seq.count("TCG"); TGC = seq.count("TGC")

    cpg1 = (2*CG + 2*GC) / (n - 1) if n > 1 else 0
    cpg2 = (ACG + AGC + CAG + CCG + CGA + CGC +
            2*CGG + CGT + CTG + GAC + GCA + 2*GCC +
            GCG + GCT + 2*GGC + GTC + TCG + TGC) / (n - 2) if n > 2 else 0
    cpg3 = (4*CAG + CCG + CGG + 4*CTG +
            4*GAC + GCC + GGC + 4*GTC)             / (n - 2) if n > 2 else 0

    return {
        'Sequence':        sequence,
        'Kmers':           kmers,
        'Label':           label,
        'region':          region,
        'A_Fraction':      A / n,
        'C_Fraction':      C / n,
        'G_Fraction':      G / n,
        'T_Fraction':      T / n,
        'PurPyr_Fraction': (A + G - C - T) / n,
        'AmKe_Fraction':   (A + C - G - T) / n,
        'WeSt_Fraction':   (A + T - C - G) / n,
        'CpG1':            cpg1,
        'CpG2':            cpg2,
        'CpG3':            cpg3,
    }


# ══════════════════════════════════════════════════════════════════════════════
# 3 – DNAshapeR MAD FEATURES  (via rpy2)
# ══════════════════════════════════════════════════════════════════════════════

def _write_fasta(sequences: list, fasta_path: str):
    """Write a minimal FASTA file."""
    with open(fasta_path, 'w') as fh:
        for i, seq in enumerate(sequences):
            fh.write(f">seq{i}\n{seq.upper()}\n")


def _diffsum(vec: np.ndarray) -> float:
    """
    Absolute cumulative difference across consecutive positions:
      sum( |pos[i+1] - pos[i]| )  for i in 0..n-2
    NaN positions are dropped before differencing.
    Matches the R rowSums(abs(t(apply(mat, 1, diff)))) logic.
    """
    v = vec[~np.isnan(vec)]
    if len(v) < 2:
        return np.nan
    return float(np.sum(np.abs(np.diff(v))))   # sum of |v[i+1] - v[i]|


def compute_shape_diffsum(sequences: list) -> pd.DataFrame:
    """
    Given a list of DNA sequences, return a DataFrame  (n_seqs x 14)
    with columns  <Feature>_diffsum  for each of the 14 DNAshapeR features.
    Falls back to NaN columns if rpy2 / DNAshapeR is unavailable.
    """
    nan_df = pd.DataFrame(
        {f"{f}_diffsum": [np.nan] * len(sequences) for f in SHAPE_FEATURES}
    )

    if not RPY2_AVAILABLE:
        print("  ⚠️  Skipping shape features (rpy2 unavailable).")
        return nan_df

    try:
        # Set personal R library path before any importr calls
        ro.r('.libPaths(c("~/R/library", .libPaths()))')

        dnashaper = importr('DNAshapeR')

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "seqs.fa")
            _write_fasta(sequences, fasta_path)

            r_shapes = dnashaper.getShape(
                fasta_path,
                shapeType=StrVector(SHAPE_FEATURES)
            )

            result = {}
            for feat in SHAPE_FEATURES:
                col = f"{feat}_diffsum"
                try:
                    mat = np.array(r_shapes.rx2(feat))
                    if mat.ndim == 1:
                        mat = mat.reshape(len(sequences), -1)
                    result[col] = [_diffsum(mat[i]) for i in range(len(sequences))]
                except Exception as e:
                    print(f"    ⚠️  Could not extract {feat}: {e}")
                    result[col] = [np.nan] * len(sequences)

        return pd.DataFrame(result)

    except Exception as e:
        print(f"  ❌ DNAshapeR call failed: {e}")
        return nan_df


# ══════════════════════════════════════════════════════════════════════════════
# 4 – WRITE SPLIT FILES
# ══════════════════════════════════════════════════════════════════════════════

def write_split(df_split: pd.DataFrame, split_name: str, outdir: str):
    """
    <split>.tsv          →  Sequence (k-mers), Label          [DNABERT input]
    sequence_<split>.tsv →  Sequence (k-mers), Label + all 24 feature columns
                             (10 composition + 14 shape diffsum; raw nucleotide excluded)

    'Kmers' is renamed to 'Sequence' in both files.
    Row i in <split>.tsv == row i in sequence_<split>.tsv.
    """
    out_df = df_split.drop(columns=['Sequence']).rename(columns={'Kmers': 'Sequence'})

    split_path = os.path.join(outdir, f"{split_name}.tsv")
    prop_path  = os.path.join(outdir, f"sequence_{split_name}.tsv")

    out_df[['Sequence', 'Label']].to_csv(split_path, sep="\t", index=False)
    out_df.to_csv(prop_path, sep="\t", index=False)

    print(f"  ✅ {split_name}.tsv           ({len(out_df)} rows, 2 cols)")
    print(f"  ✅ sequence_{split_name}.tsv  ({len(out_df)} rows, {len(out_df.columns)} cols)")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    if len(sys.argv) != 4:
        print("Usage: python pipeline.py <region> <infile.tsv> <outdir>")
        sys.exit(1)

    region = sys.argv[1]
    infile = sys.argv[2]
    outdir = sys.argv[3].rstrip("/")
    os.makedirs(outdir, exist_ok=True)

    # ── 1. Load & validate ────────────────────────────────────────────────────
    print(f"\n[1] Loading '{infile}' ...")
    raw = pd.read_csv(infile, sep="\t")

    col_map   = {c.lower(): c for c in raw.columns}
    seq_col   = col_map.get('sequence', col_map.get('seq', None))
    label_col = col_map.get('label',    col_map.get('class', col_map.get('target', None)))

    if seq_col is None or label_col is None:
        print(f"  ❌ Could not find Sequence/Label columns. Found: {list(raw.columns)}")
        sys.exit(1)

    print(f"  Using columns: sequence='{seq_col}', label='{label_col}'")
    df = raw[[seq_col, label_col]].rename(columns={seq_col: 'Sequence', label_col: 'Label'}).copy()
    df['Sequence'] = df['Sequence'].str.upper()

    before = len(df)
    df = df[df['Sequence'].str.fullmatch(r'[ACGT]+', na=False)].reset_index(drop=True)
    dropped = before - len(df)
    if dropped:
        print(f"  ⚠️  Dropped {dropped} rows with non-ACGT / missing sequences.")
    print(f"  {len(df)} valid rows  |  Label distribution:\n"
          f"{df['Label'].value_counts().to_string()}")

    # ── 2. Generate 6-mers ────────────────────────────────────────────────────
    print(f"\n[2] Generating {K}-mers (sliding window = 1) ...")
    df['Kmers'] = df['Sequence'].apply(lambda s: generate_kmers(s, K))

    # ── 3. Stratified 70 / 15 / 15 split ─────────────────────────────────────
    print(f"\n[3] Stratified split  70 / 15 / 15 ...")
    train_df, temp_df = train_test_split(
        df,
        test_size=(1 - TRAIN_RATIO),
        stratify=df['Label'],
        random_state=RANDOM_SEED,
    )
    test_df, dev_df = train_test_split(
        temp_df,
        test_size=DEV_RATIO / (TEST_RATIO + DEV_RATIO),
        stratify=temp_df['Label'],
        random_state=RANDOM_SEED,
    )
    splits = {
        'train': train_df.reset_index(drop=True),
        'test':  test_df.reset_index(drop=True),
        'dev':   dev_df.reset_index(drop=True),
    }
    for name, sdf in splits.items():
        print(f"  {name:5s}: {len(sdf):>6,} rows")

    # ── 4. Composition features ───────────────────────────────────────────────
    print(f"\n[4] Computing sequence composition features ...")
    for name, sdf in splits.items():
        records = [
            compute_composition_features(
                row['Sequence'], row['Label'], row['Kmers'], region
            )
            for _, row in sdf.iterrows()
        ]
        splits[name] = pd.DataFrame(records)

    # ── 5. DNAshapeR MAD features ─────────────────────────────────────────────
    print(f"\n[5] Computing DNAshapeR diffsum features (14 shape parameters) ...")
    for name, sdf in splits.items():
        print(f"  → {name} split ({len(sdf):,} sequences) ...")
        shape_df = compute_shape_diffsum(sdf['Sequence'].tolist())
        splits[name] = pd.concat(
            [sdf.reset_index(drop=True), shape_df.reset_index(drop=True)],
            axis=1
        )
    print(f"\n  Final feature columns ({len(splits['train'].columns)}):")
    print("  " + ", ".join(splits['train'].columns.tolist()))

    # ── 6. Write output files ─────────────────────────────────────────────────
    print(f"\n[6] Writing output files → '{outdir}/' ...")
    for name, sdf in splits.items():
        write_split(sdf, name, outdir)

    print("\n✅ Pipeline complete.\n")


if __name__ == "__main__":
    main()

# python /data/projects/dna/DNABERT/examples/scripts/props.py Promoter /data/projects/CorePromoterPWM/epd/dnabert_corepromMM_package/Promoter_nonpromoter.tsv /data/projects/CorePromoterPWM/epd/dnabert_corepromMM_package/data_props/
# python /data/projects/dna/DNABERT/examples/scripts/props.py Promoter /data/projects/CorePromoterPWM/epd/dnabert_corepromMM_package/Promoter_nonpromoter.tsv /data/projects/CorePromoterPWM/epd/dnabert_corepromMM_package/data_props/

def run(region: str, infile: str, outdir: str, r_lib: str = None):
    import sys
    sys.argv = ["coreprom-props", region, infile, outdir]
    main()
