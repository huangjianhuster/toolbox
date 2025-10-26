#!/usr/bin/env python3
"""
Simple argparse template with various argument types.
"""

import argparse


def main(args):
    """
    Main function to process command line arguments.

    Args:
        args: Parsed command line arguments
    """
    print("=== Parsed Arguments ===")
    print(f"Integer value: {args.int_val}")
    print(f"Float value: {args.float_val}")
    print(f"String value: {args.str_val}")
    print(f"Multiple values: {args.multi_val}")
    print(f"Colon-separated values: {args.colon_val}")

    # Example: Process colon-separated values
    if args.colon_val:
        print("\n=== Parsed Colon-Separated Values ===")
        for item in args.colon_val:
            parts = item.split(':')
            print(f"  {item} -> {parts}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Template script with various argparse argument types"
    )

    # Single value arguments with default types (always required)
    parser.add_argument(
        "-i", "--int-val",
        type=int,
        default=10,
        required=True,
        help="Integer value (default: 10)"
    )

    parser.add_argument(
        "-f", "--float-val",
        type=float,
        default=3.14,
        required=True,
        help="Float value (default: 3.14)"
    )

    parser.add_argument(
        "-s", "--str-val",
        type=str,
        default="example",
        required=True,
        help="String value (default: 'example')"
    )

    # Multiple value input (default: None, always required)
    parser.add_argument(
        "-m", "--multi-val",
        nargs="+",
        default=None,
        required=True,
        help="Multiple values (space-separated)"
    )

    # Multiple colon-separated values (default: None, always required)
    parser.add_argument(
        "-c", "--colon-val",
        nargs="+",
        default=None,
        required=True,
        help="Multiple colon-separated values (e.g., key1:value1 key2:value2)"
    )

    # Boolean flag argument (default: False, becomes True when provided)
    parser.add_argument(
        "-f", "--flag",
        action="store_true",
        help="Enable flag (default: False)"
    )

    args = parser.parse_args()
    main(args)
