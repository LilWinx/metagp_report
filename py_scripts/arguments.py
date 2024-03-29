import argparse

__version__ = "0.0.1"


def create_parser():
    parser = argparse.ArgumentParser(description="Meta-GP_Report")
    parser.add_argument("--wgsid", "-w", required=True, help="Sample WGS ID")
    parser.add_argument("--output", "-o", help="Output directory to write to")
    parser.add_argument("--input", "-i", help="Input folder - Autoread all required files to generate output")
    parser.add_argument("--natype", "-n", choices=[
            "dna",
            "rna",
            "both",
        ],
        required=True, 
        help="Select NA type of sample")
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        help="get Meta-GP_Report version",
        version=f"Meta-GP_Report v{__version__}")
    return parser
