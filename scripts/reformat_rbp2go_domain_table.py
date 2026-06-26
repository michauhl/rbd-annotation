#!/usr/bin/env python3

import argparse
import csv
import sys

"""
python reformat_rbp2go_domain_table.py \
    -i ../data/Table_selected_new_RBDs.txt \
    -o ../data/rbp2go.selected_new_rbds.reform.txt
    
"""

def find_header_line(infile, expected_first_col="InterPro_ID"):
    """
    Skip metadata lines until the real quoted TSV header is found.
    Returns a csv reader positioned after the header.
    """
    for line_num, line in enumerate(infile, start=1):
        stripped = line.strip()

        if not stripped:
            continue

        fields = next(csv.reader([stripped], delimiter="\t", quotechar='"'))

        if fields and fields[0] == expected_first_col:
            return fields, line_num

    sys.exit(f"Error: could not find header line starting with '{expected_first_col}'.")


def main():
    parser = argparse.ArgumentParser(
        description="Convert quoted RBP2GO-style InterPro table to clean TSV."
    )
    parser.add_argument("-i", "--input", required=True, help="Input table file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    with open(args.input, "r", newline="") as infile, \
         open(args.output, "w", newline="") as outfile:

        header, header_line_num = find_header_line(infile)

        reader = csv.reader(infile, delimiter="\t", quotechar='"')
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")

        writer.writerow(header)

        n_rows = 0

        for line_num, row in enumerate(reader, start=header_line_num + 1):
            if not row or all(cell.strip() == "" for cell in row):
                continue

            if len(row) != len(header):
                sys.exit(
                    f"Error: line {line_num} has {len(row)} columns, "
                    f"expected {len(header)}."
                )

            writer.writerow(row)
            n_rows += 1

    print(f"Wrote {n_rows} rows to: {args.output}")


if __name__ == "__main__":
    main()
