#!/usr/bin/env python3
"""
Generate an IGV session based on output of a pb-human-wgs-workflow snakemake.
"""

__author__ = "William Rowell"
__version__ = "0.1.0"


import jinja2
import argparse


def main(args):
    ## Render template
    env = jinja2.Environment()
    env.loader = jinja2.FileSystemLoader(".")
    template = env.get_template("templates/cohortigvxml.jinja")
    figv = open(args.output, "w")
    figv.write(template.render(
        webroot=args.webroot,
        reference=args.reference,
        cohort=args.cohort,
        samples=args.samples.split(","),
        segdups=args.segdups,
        oddregions=args.oddregions,
        repeats=args.repeats
        ))
    figv.close()


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("webroot", help="Root URL.")
    parser.add_argument("reference", help="Short reference name used in filenames.")
    parser.add_argument("cohort", help="Cohort ID")
    parser.add_argument("samples", help="Comma delimited list of samples.")
    parser.add_argument("segdups", help="Relative path to SegDup bed.")
    parser.add_argument("oddregions", help="Relative path to OddRegions bed.")
    parser.add_argument("repeats", help="Relative path to Repeats bed.")
    parser.add_argument("output", help="IGV session xml.")
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)