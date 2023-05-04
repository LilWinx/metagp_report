import argparse

__version__ = "0.0.1"


def create_parser():
    parser = argparse.ArgumentParser(description="Meta-GP_Report")
    parser.add_argument("--outdir", "-o", help="Output directory to write to")
    parser.add_argument("--txt", "-t", help="Text file input", required=True)
    parser.add_argument("--krona", "-k", help="Path to R2 file", required=True)
    parser.add_argument("--manual", "-m", help="Enter patient details manually",)
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        help="get Meta-GP_Report version",
        version=f"Meta-GP_Report v{__version__}")
    return parser
