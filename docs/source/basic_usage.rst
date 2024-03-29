Basic Usage
===========

`Example 1`_ (below) illustrates usage of the ADCIRCpy Python API to generate a "best track" parametric wind configuration for running ADCIRC on hurricane Sandy (2012).

.. _Example 1:

.. code-block:: python

    #!/usr/bin/env python

    from datetime import timedelta

    from adcircpy import AdcircMesh, AdcircRun, Tides
    from adcircpy.forcing.winds import BestTrackForcing
    from adcircpy.server import SlurmConfig

    # load an ADCIRC mesh grid from a `fort.14` file to a new mesh object
    mesh = AdcircMesh.open('fort.14', crs='epsg:4326')

    # add nodal attributes from a `fort.13` file to the mesh object
    mesh.import_nodal_attributes('fort.13')

    # create a tidal forcing object, using all constituents
    tidal_forcing = Tides()
    tidal_forcing.use_all()

    # add data from the tidal forcing object to the mesh object
    mesh.add_forcing(tidal_forcing)

    # create a wind forcing object for Hurricane Sandy (2012)
    wind_forcing = BestTrackForcing('Sandy2012')

    # add wind forcing data to the mesh object
    mesh.add_forcing(wind_forcing)

    # create a Slurm (HPC job manager) configuration object.
    slurm = SlurmConfig(
        account='account',
        ntasks=1000,
        run_name='ADCIRCpy documentation example',
        partition='partition',
        walltime=timedelta(hours=8),
        mail_type='all',
        mail_user='example@email.gov',
        log_filename='example.log',
        modules=['intel/2020', 'impi/2020', 'netcdf/4.7.2-parallel'],
        path_prefix='$HOME/adcirc/build',
    )

    # create an ADCIRC run driver object
    driver = AdcircRun(
        mesh=mesh,
        server_config=slurm,
        spinup_time=timedelta(days=15),
    )

    # write configuration files to the specified directory
    driver.write(output_directory="./model_inputs")

The set of files generated by this configuration (``fort.14``, ``fort.15.coldstart``, ``fort.15.hotstart``, ``fort.22``, and ``slurm.job``) represents the minimum working example for an ADCIRC run using parametric winds. To submit this configuration to the Slurm job manager, run the following command in a shell with access to the ADCIRC binaries:

.. code-block:: shell

    sbatch slurm.job

Note that this setup generates model inputs based on default values. These may not be necessarily optimal. For this reason, it is advisable that users perform checks of the outputs and introduce optimizations by modifying the parameters given to the Python API. For example, the user might want to check that the auto-computed timestep matches the expectations from experience. Instead of doing manual modifications to the ``fort.15`` file, the user makes use of the ADCIRCpy module to customize and optimize their modeling needs. If the user, for example, requires their model to use a specific timestep of two seconds, it would suffice to execute the following statement:

.. code-block:: python

    driver.timestep = 2

before :code:`driver.write()`. Most of the ``fort.15`` options can be overridden directly by the user, with a very small exception of parameters that should be considered "private" in the context of the ``fort.15`` file.

.. code-block:: python

    # Modify timestep to 2 seconds.
    driver.timestep = 2.

    # Add a constant Mannings coefficient field to the mesh
    mesh.mannings_n_at_sea_floor = mesh.coords.shape[0]*[0.025]

    # generate TAU0 factors
    mesh.generate_tau0()

    # Write new model configuration to disk.
    driver.write("model_inputs_modified", overwrite=True)

After these modifications, the resulting directory contains ``fort.13``, ``fort.14``, ``fort.15.coldstart``, ``fort.15.hotstart``, ``fort.22``, and ``slurm.job``. The new ``fort.13`` includes the newly-added Manning's N and :code:`TAU0` factors.
