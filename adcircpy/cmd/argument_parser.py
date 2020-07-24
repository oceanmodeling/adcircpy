import argparse
import logging
import sys


def mesh(parser):

    # mesh
    parser.add_argument('mesh')

    # Mesh spatial reference
    msg = "Mesh projection information. Defaults to EPSG:4326 which "
    msg += "corresponds to WGS84. For Cartesian meshes use 3395 which "
    msg += "corresponds to Mercator projection."
    parser.add_argument('--crs')


def output_directory(parser):
    # output directory
    msg = "Directory to which ADCIRC input files will be written to. "
    parser.add_argument('--output-directory', '--outdir', "-o", help=msg)


def allow_overwrite(parser):
    msg = 'Allows overwrite of output directory.'
    parser.add_argument('--overwrite', help=msg, action="store_true")


def generate_only(parser):
    msg = "Generates and saves input files to the output directory but does "
    msg += "not deploy the ADCIRC run."
    parser.add_argument(
        '--generate-only', "--no-run", "--skip-run",
        action="store_true",
        help=msg
        )


def log_level(parser):
    log_level = parser.add_mutually_exclusive_group()
    log_level.add_argument(
        '--log-level-info',
        nargs='?',
        const=logging.INFO,
        dest="log_level")
    log_level.add_argument(
        '--log-level-debug',
        nargs='?',
        const=logging.DEBUG,
        dest="log_level")
    log_level.add_argument(
        '--log-level-warning',
        nargs='?',
        const=logging.WARNING,
        dest="log_level")


def server(parser):

    # flag some options as required when a resource manager is enabled
    _required = "--use-torque" in sys.argv
    _required = _required | ("--use-pbs" in sys.argv)
    _required = _required | ("--use-slurm" in sys.argv)

    # add server options
    parser.add_argument('--hostname')
    parser.add_argument('--port', type=int)
    parser.add_argument("--wdir", required=_required)
    parser.add_argument("--keep-wdir", action="store_true")
    parser.add_argument(
        "--binaries-path", "--binaries-prefix",
        dest="binaries_prefix")
    parser.add_argument("--source-script")
    parser.add_argument("--additional-mpi-options")

    # make nproc required when using ssh
    args = parser.parse_known_args()[0]
    if args.hostname is not None:
        parser.add_argument("--nproc", "--ncpu", type=int, required=True)
    else:
        parser.add_argument("--nproc", "--ncpu", type=int, default=-1)

    # add resource manager option
    manager = parser.add_mutually_exclusive_group()
    manager.add_argument('--use-torque', action="store_true")
    manager.add_argument('--use-pbs', action="store_true")
    manager.add_argument('--use-slurm', action="store_true")

    # resource manager specific options
    parser.add_argument('--account', required=_required)
    parser.add_argument('--walltime', required=_required)
    parser.add_argument(
        '--module',
        default=list(),
        action='append',
        dest='modules'
        )


def tidal_constituents(parser):
    # tidal constituents
    msg = "Tidal constituent to be forced in the model. Pass "
    msg += "--use-constituent='all' to use all available constituents "
    msg += "(K1, O1, P1, Q1, MM, Mf, M4, MN4, MS4, 2N2, S1) "
    msg += "or use --use-constituent='major'. "
    msg += "For a custom list of forcing constituents, pass -c= for each "
    msg += "individual constituent to use (case-insensitive). "
    msg += "Use None for no tidal forcing. Defaults to 'all'."
    parser.add_argument(
        "--constituents", "-c",
        action='append',
        nargs="?",
        choices=["K1", "O1", "P1", "Q1", "MM", "Mf", "M4", "MN4", "MS4",
                 "2N2", "S1", "all", "major"],
        dest='constituents',
        default=[],
        help=msg
        )


def timestep(parser):
    parser.add_argument("--timestep", type=float, default=False)


def gwce_solution_scheme(parser):
    parser.add_argument(
        "--gwce-solution-scheme",
        choices=['semi-implicit', 'explicit'],
        # default='explicit'
        default='semi-implicit'
        )


def boundaries_generation(parser):
    parser.add_argument("--generate-boundaries", action="store_true")
    parser.add_argument("--boundaries-threshold", default=0., type=float)
    parser.add_argument(
        "--land-ibtype",
        default=20,
        type=int,
        choices=[0, 10, 20]
        )
    parser.add_argument(
        "--island-ibtype",
        default=21,
        type=int,
        choices=[1, 11, 21]
        )


def best_track(parser):
    # storm_id
    msg = "National Hurricane Center (NHC) storm id. "
    msg += " Examples: AL132012 for Sandy2012 or AL152017 for Maria2017."
    parser.add_argument('storm_id', help=msg)
    parser.add_argument('--start-date')
    parser.add_argument('--end-date')
    parser.add_argument('--spinup-days', type=float, required=True)


def tidal_run(parser):
    # start_date
    msg = "Start date is relative to hotstart, that is, this is the "
    msg += "true start date of the model (in UTC time). Use format "
    msg += "%%Y-%%m-%%dT%%H:%%M to specify date. For example, for August "
    msg += "1, 2013, 00:00 hours, write \"2018-08-01T00:00\" (can be used "
    msg += "with or without the quotes)."
    parser.add_argument('start_date', help=msg)
    # end_date
    parser.add_argument('end_date')
    # spinup_days
    parser.add_argument('--spinup-days', type=float, required=True)


def timezone(parser):
    parser.add_argument("--timezone")


def netcdf(parser):
    parser.add_argument(
        '--ascii',
        dest='netcdf',
        action='store_false',
        default=True,
        help="Request outputs in ASCII format. NetCDF is the default."
    )


def tau0(parser):
    # TAU0 override
    parser.add_argument("--generate-tau0", "--tau0-gen", action="store_true")

    # FFACTOR override
    msg = "This is the constant friction coefficient used if any of "
    msg += "'quadratic_friction_coefficient_at_sea_floor', "
    msg += "'mannings_n_at_sea_floor', "
    msg += "'chezy_friction_coefficient_at_sea_floor', "
    msg += "'bottom_roughness_length', "
    msg += "are not passed as a nodal attribute. "
    msg += "Defaults to 0.02."
    parser.add_argument("--FFACTOR", type=float, default=0.02, help=msg)


def nodal_attributes(parser):
    parser.add_argument("--fort13")
    # attributes
    msg = "Use coldstart attribute that exists in fort.13 during coldstart "
    msg += "phase. Use --coldstart-attribute=all to use all available "
    msg += "attributes."
    parser.add_argument(
        "--coldstart-attribute",
        default=['all'],
        action='append',
        dest='coldstart_attributes',
        help=msg
        )

    # hotstart attributes
    msg = "Use nodal attribute that exists in fort.13 during hotstart "
    msg += "phase. Use all to use all available attributes or use None."
    parser.add_argument(
        "--hotstart-attribute",
        default=['all'],
        action='append',
        dest='hotstart_attributes',
        help=msg)


def surface_output(physical_var, parser, spinup=False):
    # surface output requests
    long_name = f"--{physical_var}-surface-sampling-rate"
    short_name = f"--{physical_var[:4]}"
    if spinup:
        long_name += "-spinup"
        short_name += "-spinup"
    msg = f"{physical_var.capitalize()} surface output sampling frequency "
    msg += f"(in minutes). When this number is greater than 0, {physical_var} "
    msg += " surface outputs are written to disk during hotstart phase."
    parser.add_argument(
        long_name,
        short_name,
        type=float,
        help=msg
        )

    # # surface output start
    # msg = f"Set {physical_var} surface output starting time in days (after "
    # msg += "model start_time)"
    # parser.add_argument(
    #     f"--{physical_var}-surface-start",
    #     f"--{physical_var[:4]}-start",
    #     type=float,
    #     help=msg
    #     )

    # # surface output end
    # msg = ""
    # parser.add_argument(
    #     f"--{physical_var}-surface-end",
    #     f"--{physical_var[:4]}-end",
    #     type=float,
    #     help=msg
    #     )
    if not spinup:  # disable harm for spinup, less clutter.
        long_name = f'--{physical_var}-surface-harmonic-analysis'
        short_name = f"--{physical_var[:4]}-harm"
        if spinup:
            long_name += "-spinup"
            short_name += "-spinup"
        # surface harmonic analysis
        msg = f"Enables {physical_var} surface harmonic analysis."
        parser.add_argument(
            long_name,
            short_name,
            action='store_true',
            default=False,
            help=msg
            )


def stations_output(physical_var, parser, spinup=False):
    long_name = f"--{physical_var}-stations-sampling-rate"
    short_name = f"--{physical_var[:4]}-stat"
    if spinup:
        long_name += "-spinup"
        short_name += "-spinup"
    # stations output requests
    msg = f"{physical_var.capitalize()} stations sampling frequency in "
    msg += f"minutes. When this number is greater than 0, {physical_var} "
    msg += "stations output is turned on during hotstart phase."
    parser.add_argument(
        long_name,
        short_name,
        type=float,
        help=msg
    )

    # # stations start
    # msg = f"Set {physical_var} stations output starting time in days (after "
    # msg += "model start_time)"
    # parser.add_argument(
    #     f"--{physical_var}-stations-start",
    #     f'--{physical_var[:4]}-s-start',
    #     type=float,
    #     help=msg
    #     )

    # # stations end
    # msg = f"Set {physical_var} stations output end time in days (after model "
    # msg += "start_time)"
    # parser.add_argument(
    #     f"--{physical_var}-stations-end",
    #     f'--{physical_var[:4]}-s-end',
    #     type=float,
    #     help=msg
    #     )
    if not spinup:  # disable harm for spinup, less clutter.
        long_name = f'--{physical_var}-stations-harmonic-analysis'
        short_name = f"--{physical_var[:4]}-s-harm"
        if spinup:
            long_name += "-spinup"
            short_name += "-spinup"
        # stations harmonic analysis
        msg = f"Enables {physical_var} stations harmonic analysis."
        parser.add_argument(
            long_name,
            short_name,
            action='store_true',
            default=False,
            help=msg
            )


def spinup_outputs(parser):
    surface_output('elevation', parser, spinup=True)
    surface_output('velocity', parser, spinup=True)
    surface_output('meteorological', parser, spinup=True)
    surface_output('concentration', parser, spinup=True)
    stations_output('elevation', parser, spinup=True)
    stations_output('velocity', parser, spinup=True)
    stations_output('meteorological', parser, spinup=True)
    stations_output('concentration', parser, spinup=True)


def outputs(parser):

    # add surface output requests
    surface_output('elevation', parser)
    surface_output('velocity', parser)
    surface_output('meteorological', parser)
    surface_output('concentration', parser)

    # parse stations from a file
    msg = "File containing list of stations for outputs. It will parse "
    msg += "the stations below the NOUTE, NOUTV and NOUTM keywords for "
    msg += "their respective stations list."
    parser.add_argument(
        "--stations-file",
        help=msg
        )

    # stations output requests
    stations_output('elevation', parser)
    stations_output('velocity', parser)
    stations_output('meteorological', parser)
    stations_output('concentration', parser)

    # spinup options
    spinup_outputs(parser)


def get_parser(runtype=None, description=None):
    parser = argparse.ArgumentParser(description=description)
    mesh(parser)
    if runtype is not None:
        if runtype == "tidal":
            tidal_run(parser)
        elif runtype == 'best_track':
            best_track(parser)
    nodal_attributes(parser)
    timezone(parser)
    tau0(parser)
    tidal_constituents(parser)
    timestep(parser)
    gwce_solution_scheme(parser)
    boundaries_generation(parser)
    output_directory(parser)
    allow_overwrite(parser)
    generate_only(parser)
    outputs(parser)
    netcdf(parser)
    log_level(parser)
    server(parser)
    return parser
