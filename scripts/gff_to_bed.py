#!/usr/bin/env python3
import argparse
import sys


def parse_attrs(attr_str: str) -> dict:
    attrs = {}
    for field in attr_str.strip().split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        elif " " in field:
            k, v = field.split(" ", 1)
        else:
            k, v = field, ""
        attrs[k.strip()] = v.strip().strip('"')
    return attrs


def main():
    ap = argparse.ArgumentParser(
        description="Convert GFF3/GTF-like annotations to BED (BED3+1: chrom, start, end, name)."
    )
    ap.add_argument("gff", help="Input GFF3/GTF file (must match your reference)")
    ap.add_argument("out_bed", help="Output BED file path")
    ap.add_argument(
        "--feature",
        default="gene",
        help="Feature type to include (GFF column 3), default: gene",
    )
    ap.add_argument(
        "--name-attr",
        default="Name",
        help="GFF attribute key used as the BED name (e.g., Name, gene, ID). Default: Name",
    )
    ap.add_argument(
        "--names",
        help="Optional file with line-delimited names to include (matches chosen name-attr)",
    )
    ap.add_argument(
        "--regex",
        help="Optional regex to include names (applied to chosen name-attr)",
    )
    args = ap.parse_args()

    include_set = None
    if args.names:
        include_set = set(
            x.strip() for x in open(args.names, "r", encoding="utf-8") if x.strip()
        )

    name_re = None
    if args.regex:
        import re

        name_re = re.compile(args.regex)

    written = 0
    with open(args.gff, "r", encoding="utf-8", errors="replace") as fin, open(
        args.out_bed, "w", encoding="utf-8"
    ) as fout:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts
            if args.feature and ftype != args.feature:
                continue
            try:
                start0 = max(0, int(start) - 1)  # BED is 0-based
                end1 = int(end)  # half-open
            except ValueError:
                continue
            name = ""
            a = parse_attrs(attrs)
            if args.name_attr in a:
                name = a[args.name_attr]
            elif "Name" in a:
                name = a["Name"]
            elif "ID" in a:
                name = a["ID"]

            if include_set is not None and name not in include_set:
                continue
            if name_re is not None and not name_re.search(name):
                continue

            fout.write(f"{chrom}\t{start0}\t{end1}\t{name}\n")
            written += 1

    if written == 0:
        print("WARNING: no features written to BED (check filters)", file=sys.stderr)


if __name__ == "__main__":
    main()

