#! /usr/bin/env python
import argparse
from adcircpy import AdcircMesh


def main():
    args = parse_args()

    # read mesh
    mesh = AdcircMesh.open(args.mesh)
    mesh.triplot(show=True)


def parse_args():
    # init arparse
    parser = argparse.ArgumentParser(
        description="Program to generate ADCIRC input files for a given mesh "
        + "hindcast storm name.")

    # mesh
    parser.add_argument(
        "mesh",
        help="Path to the mesh file to be used.")

    # storm id
    parser.add_argument(
        "storm_id",
        help="National Hurricane Center storm id. Example: AL132012 for "
        + "Sandy 2012 or AL152017 for Maria 2017.")

    # output path
    parser.add_argument(
        "output_path",
        help="Directory where to save output files. If the directory does not "
        + "exist it will be created.")

    # overwrite flag
    parser.add_argument(
        "--overwrite",
        help="Optional flag to allow for output file overwrite.")

    # parser.add_argument(
    #     '--start-date', help="format is %Y%m%d%H")

    # parser.add_argument(
    #     '--end-date', help="format is %Y%m%d%H")

    return parser.parse_args()


if __name__ == "__main__":
    main()
