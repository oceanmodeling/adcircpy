#! /usr/bin/env python
import argparse
from datetime import datetime, timedelta
import pathlib

from adcircpy import AdcircMesh, Fort15
from adcircpy.forcing.tides import Tides


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('mesh')
    parser.add_argument('start_date', type=lambda x: datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))
    parser.add_argument('run_days', type=float)
    parser.add_argument('--output-file', type=pathlib.Path)
    parser.add_argument('--tidal-database', '--tidal-db', choices=['hamtide', 'tpxo'])
    parser.add_argument('--mesh-crs')
    return parser.parse_args()


def main():
    args = parse_args()
    tides = Tides(tidal_source=args.tidal_database)
    tides.use_all()
    tides.start_date = args.start_date
    tides.end_date = tides.start_date + timedelta(days=args.run_days)
    mesh = AdcircMesh.open(args.mesh, crs=args.mesh_crs)
    mesh.add_forcing(tides)
    fort15 = Fort15(mesh)
    if args.output_file is not None:
        with open(args.output_file, 'w') as f:
            f.write(fort15.get_tidal_forcing())
    else:
        print(fort15.get_tidal_forcing())


if __name__ == '__main__':
    main()
