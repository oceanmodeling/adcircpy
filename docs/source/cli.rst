CLI Commands
============

The :code:`adcircpy` command line interface (CLI) provides commands accessible from the command line, allowing a smart way to produce forecasts, tidal-only hindcasts, "best track" hindcast operations, and practicing parametric winds appropriate for storm surge hindcasts, for rapid iteration and debugging. A single command generates a sanitized file through Python code, which is user-customized to the virtual model scenario.

``tide_gen``
------------
The :code:`tide_gen` command generates the tidal forcing table required in the ``fort.15`` file. The generated table is suitable for being pasted onto a target ``fort.15`` file. In reality, this algorithm calls the same functions invoked when generating full model runs, so it can be used to check whether the tidal forcing outputs will look reasonable before running the entirety of the algorithm.

.. code-block:: shell

    tide_gen /path/to/your/fort.14 '2021-02-26T00:00:00' 15 --mesh-crs='epsg:4326'

.. program-output:: tide_gen -h

``tidal_run``
-------------
The :code:`tidal_run` entrypoint will generate the necessary set of input files for forecast and hindcast tidal-only runs with any mesh. Most of the options that tidal_run takes are also valid for other command line entry points, so getting familiar with the options on tidal_run is helpful to understand the options for more advanced CLI entry points. A tidal only run is perhaps the first step towards mesh validation, and mesh stability checks, therefore having an entrypoint that is able to generate ADCIRC input files for tidal only runs pays greatly for initial mesh stability checks.

.. code-block:: shell

    tidal_run \
        /path/to/fort.14 \
        $(date +"%Y-%m-%dT%H:%M:%S") \
        15 \
        --spinup-days=5  \
        --tau0-gen  \
        --crs=EPSG:4326  \
        --constituents=all \
        --timestep=10. \
        --stations-file=stations.txt \
        --elev-stat=6. \
        --generate-linear-mannings

.. program-output:: tidal_run -h

``best_track_run``
------------------
The :code:`best_track_run` command line entry point was designed to be able to generate a full parametric wind run in ADCIRC from the command line. This entry point essentially encapsulates most of the :code:`adcircpy` functionality in a single command. The following is an example of how to invoke this functionality from the bash shell.

.. code-block:: shell

    best_track_run \
        fort.14 \
        Ike2008 \
        --spinup-days=15 \
        --crs=EPSG:4326 \
        --fort13=fort.13 \
        --output-directory=Ike2008 \
        --constituents=major \
        --skip-run \
        --generate-linear-mannings \
        --tau0-gen \
        --timestep=1.0 \
        --elev=60. \
        --stations-file=coops.txt \
        --elev-stat=6. \
        --binaries-prefix=.python_env/bin \
        --overwrite \
        --use-slurm \
        --account=nosofs \
        --slurm-ntasks=800 \
        --partition=orion \
        --walltime=8 \
        --mail-type=all \
        --mail-user=jaime.calzada@noaa.gov \
        --module=intel/2020 \
        --module=impi/2020 \
        --module=netcdf/4.7.2-parallel \
        --log-level=info

.. program-output:: best_track_run -h

``best_track_file``
-------------------
The :code:`best_track_file` entry point generates an ``aswip``-ready "best track" file. This uses the :code:`adcircpy.forcing.winds.BestTrackForcing` class.

.. code-block:: shell

    best_track_file Sandy2012

.. program-output:: best_track_file -h

``fort63``
----------
.. program-output:: fort63 -h

``plot_maxele``
---------------
.. program-output:: plot_maxele -h

``plot_fort61``
---------------
.. code-block:: shell

    plot_fort61 /path/to/fort.61.nc MSL --show --coops-only

.. program-output:: plot_fort61 -h

``plot_mesh``
-------------
.. code-block:: shell

    plot_mesh /path/to/fort.14 --show-elements

.. program-output:: plot_mesh -h
