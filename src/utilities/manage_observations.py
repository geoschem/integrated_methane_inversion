#!/usr/bin/env python3
"""Product-neutral commands for managing IMI observation files."""

from __future__ import annotations

import argparse
from pathlib import Path

from src.inversion_scripts.satellite_products import get_satellite_product


def parse_args() -> argparse.Namespace:
    """Parse the observation-management command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    download = subparsers.add_parser("download", help="download raw observations")
    download.add_argument("satellite_product")
    download.add_argument("start_date")
    download.add_argument("end_date")
    download.add_argument("destination", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the selected observation-management command."""
    args = parse_args()
    product = get_satellite_product(args.satellite_product)
    if args.command == "download":
        product.download(args.start_date, args.end_date, args.destination)


if __name__ == "__main__":
    main()
