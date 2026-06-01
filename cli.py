"""
CLI entry point for coreprom-MM.
Usage: coreprom-props <region> <infile.tsv> <outdir> [--r-lib PATH]
"""
import argparse
from coreprom_MM.pipeline import run


def main():
    parser = argparse.ArgumentParser(
        prog="coreprom-props",
        description=(
            "Generate 6-mer sequence files and DNA shape/composition "
            "property files from a labeled promoter TSV."
        ),
    )
    parser.add_argument("region",  help="Region label, e.g. Promoter or Enhancer")
    parser.add_argument("infile",  help="Input TSV with columns: Sequence, Label")
    parser.add_argument("outdir",  help="Output directory (created if missing)")
    parser.add_argument(
        "--r-lib", default=None, metavar="PATH",
        help="Custom R library path (default: ~/R/library)"
    )
    args = parser.parse_args()
    run(region=args.region, infile=args.infile,
        outdir=args.outdir, r_lib=args.r_lib)


if __name__ == "__main__":
    main()
