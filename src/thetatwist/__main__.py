from __future__ import annotations

import argparse
import json
from pathlib import Path

from .core import certificate, write_results


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the exact theta-twist audit")
    parser.add_argument("--d", type=int, default=3, help="odd integer d >= 3")
    parser.add_argument(
        "--write-results",
        type=Path,
        help="write the certificate and SHA-256 manifest to this directory",
    )
    args = parser.parse_args()
    result = certificate(args.d)
    print(json.dumps(result, indent=2, sort_keys=True))
    if args.write_results is not None:
        write_results(args.write_results, args.d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
