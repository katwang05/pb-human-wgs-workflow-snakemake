#!/usr/bin/env python3
"""
This script takes a vcf file and filters variants that are less likely to be causal
based on the inheritance.
"""

__author__ = "William Rowell"
__version__ = "0.1.0"


import sys
import argparse
from pysam import VariantFile


def unaffected_homalt(homalts, unaffecteds):
    "if unaff is homalt, variant is not causal"
    return homalts.intersection(unaffecteds)


def affected_homref_or_unknown(alts, affecteds):
    "if aff is homref or unk, variant is not causal"
    return alts.intersection(affecteds) != affecteds


def unaffected_hetalt_affected_hetalt(homalts, hetalts, unaffecteds, affecteds):
    "if unaff is hetalt, variant is causal only if aff is homalt"
    return (
        hetalts.intersection(unaffecteds)
        and homalts.intersection(affecteds) != affecteds
    )


def below_min_depth(record, min_depth):
    "check if any samples are below minimum depth"
    return any([record.samples[_]["DP"] < min_depth for _ in record.samples.keys()])


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("vcf_in", help="input vcf", type=str)
    parser.add_argument("vcf_out", help="output vcf", type=str)
    parser.add_argument("--affecteds", help="reference contig", type=str, default="")
    parser.add_argument("--unaffecteds", help="reference contig", type=str, default="")
    parser.add_argument("--min_depth", help="minimum GT/DP", type=int, default=0)
    args = parser.parse_args(arguments)

    affecteds = set(args.affecteds.split(",")) if args.affecteds else set()
    unaffecteds = set(args.unaffecteds.split(",")) if args.unaffecteds else set()

    with VariantFile(args.vcf_in) as reader, VariantFile(
        args.vcf_out, "w", header=reader.header
    ) as writer:
        for record in reader:
            homalts = set(record.info["homalt"]) if "homalt" in record.info else set()
            hetalts = set(record.info["hetalt"]) if "hetalt" in record.info else set()
            alts = homalts.union(hetalts)

            if unaffected_homalt(homalts, unaffecteds):
                continue
            if affected_homref_or_unknown(alts, affecteds):
                continue
            if unaffected_hetalt_affected_hetalt(
                homalts, hetalts, unaffecteds, affecteds
            ):
                continue
            if below_min_depth(record, args.min_depth):
                continue
            writer.write(record)


if __name__ == "__main__":
    """This is executed when run from the command line"""
    sys.exit(main(sys.argv[1:]))
