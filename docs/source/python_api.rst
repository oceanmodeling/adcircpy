Python API
==========

The main idea of ADCIRCpy is to encompass all necessary implementation into an object-oriented framework reusable over many different modeling scenarios. The core library is designed such that it can operate over an arbitrary mesh and avoid "hard-coding" of paths and project components. Instead, all user-defined data sets are given as inputs to relevant API components. The current discussion will be limited to the most common library components used by users. Most "private" library components will not be discussed in this document.

``adcircpy.mesh`` module
------------------------

The :code:`adcircpy.mesh` module contains all methods and functions related to input and output of an arbitrary AdcircMesh UML diagram (``fort.13``, ``fort.14``). The main user interface of the :code:`mesh` module is the :code:`adcircpy.mesh.AdcircMesh` class, which inherits from :code:`adcircpy.mesh.Fort14`. An :code:`AdcircMesh` object can read an AdcircMesh UML diagram (``fort.14``), configure model forcings, hold nodal attribute information (``fort.13``), and generate select parameters via additional methods. For example, the ``primitive_weighting_in_continuity_equation``  parameter can be generated via the :code:`.generate_tau0()` method. This method is a pythonic implementation of the :code:`tau0_gen.f90` :cite:`Weiver2008`. Another nodal attribute that can be set using a built-in method is the ``mannings_n_at_sea_floor``, which can be set via the :code:`.generate_constant_mannings_n()` or :code:`.generate_linear_mannings_n()` methods. This elucidates the convenience of the object-oriented implementation; all necessary tools exist in the same workspace and are accessible via the same object on a per-instance basis. The :code:`mesh` module also contains parsers that can read and convert between ``fort.14`` and ``*.2dm`` files.

``Grd`` class
^^^^^^^^^^^^^
The :code:`Grd` class is used as a base class that containing multiple useful functions and methods related to the specific format for an unstructured mesh in 2-dimensional space used by ADCIRC. In the context of ADCIRCpy, the :code:`Grd` class can be instantiated with a properly constructed :code:`nodes` and :code:`elements` dictionary, but it is also used as a base class for other ADCIRC objects that can be represented as conforming to a euclidean graph in 2-dimensional space (for example, the output file ``maxele.63.nc``).

.. autoclass:: adcircpy.mesh.base.Grd

`Example 2`_ shows a demonstration of how this class can be instantiated with a test graph structure, while the `Grd UML diagram`_ (`Unified Modeling Language <https://en.wikipedia.org/wiki/Unified_Modeling_Language>`_) shows the output plot of the graph.

.. _Example 2:

.. code-block:: python

    from adcircpy.mesh.base import Grd

    nodes = {
        '1': ((0., 0.), -5.),
        '2': ((.5, 0.), -4.),
        '3': ((1., 0.), -3.),
        '4': ((1., 1.), -2.),
        '5': ((0., 1.), -1.),
        '6': ((.5, 1.5), 0.),
        '7': ((.33, .33), 1.),
        '8': ((.66, .33), 2.),
        '9': ((.5, .66), 3.),
        '10': ((-1., 1.), 4.),
        '11': ((-1., 0.), 5.),
    }

    elements = {
        '1': ['5', '7', '9'],
        '2': ['1', '2', '7'],
        '3': ['2', '3', '8'],
        '4': ['8', '7', '2'],
        '5': ['3', '4', '8'],
        '6': ['4', '9', '8'],
        '7': ['4', '6', '5'],
        '8': ['5', '10', '11', '1'],
        '9': ['9', '4', '5'],
        '10': ['5', '1', '7']
    }

    fort14 = Grd(nodes, elements)

    fort14.triplot(linewidth=1, show=True)

.. figure:: figures/fort14_triplot_example.png
  :alt: Plotting result of above code

  Plotting result of `Example 2`_

.. _Grd UML diagram:

.. figure:: figures/Grd_UML.png
  :alt: UML component diagram for the :code:`Grd` class

  UML component diagram for the :code:`Grd` class

``Fort14`` class
^^^^^^^^^^^^^^^^
The :code:`Fort14` class is a subclass of the :code:`Grd` class, which is used specifically for ``fort.14`` file types.

.. autoclass:: adcircpy.mesh.fort14.Fort14

As can be seen from `Fort14 UML diagram`_ diagram, the :code:`Fort14` class essentially extends the :code:`Grd` class by adding model boundaries. While the :code:`Hull` instances in the :code:`Grd` provides geometrical boundaries for the :code:`Grd` types (which are derived from the connectivity table), the :code:`Fort14Boundaries` instance provides an API to interact with the model's boundary types: ocean, land, interior, inflow, outflow, weir and culvert.

.. _Fort14 UML diagram:

.. figure:: figures/Fort14_UML.png
  :alt: UML component diagram for the :code:`Fort14` class

  UML component diagram for the :code:`Fort14` class

``AdcircMesh`` class
^^^^^^^^^^^^^^^^^^^^
The :code:`AdcircMesh` class is the center point of the ADCIRCpy interface. It consists of a mixin class that provides an API for reading and writing AdcircMesh UML diagram files, establishing the boundary and surface forcings, as well as establishing the model's nodal attributes. The :code:`AdcircMesh` class inherits properties and methods from :code:`Fort14`, and acts as the main interface between the user and the ADCIRC model forcings and nodal attributes.

.. autoclass:: adcircpy.mesh.mesh.AdcircMesh

Experienced ADCIRC users may also recognize familiar ADCIRC attribute and method names from `AdcircMesh UML diagram`_, along with the attributes and methods inherited from :code:`Fort14`. `HSOFS mesh`_ shows the output of the :code:`AdcircMesh.make_plot()` method using the Hurricane Surge On-Demand Forecast System (HSOFS) mesh :cite:`Riverside2015`.

.. _AdcircMesh UML diagram:

.. figure:: figures/classes_adcircpy.mesh.AdcircMesh.png
  :alt: UML component diagram for the :code:`AdcircMesh` class

  UML component diagram for the :code:`AdcircMesh` class

.. _HSOFS mesh:

.. figure:: figures/hsofs_mesh.png
  :alt: HSOFS mesh plot using the :code:`AdcircMesh.make_plot()` method

  HSOFS mesh plot using the :code:`AdcircMesh.make_plot()` method

managing nodal attributes using the ``AdcircMesh`` class
""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The `AdcircMesh UML diagram`_ diagram shows that many attributes of the attributes defined in :code:`AdcircMesh` class match the nodal attribute names allowed in the ADCIRC nodal attributes file (fort.13). These attributes have been defined explicitly in the :code:`AdcircMesh` class so that the user can directly assign values to them. By doing so, these values are automatically enabled for the model and taken into account when generating the fort.15 file. The user just needs to generate the corresponding arrays. Alternatively, some nodal parameters can be populated through provided methods in the :code:`AdcircMesh` class, for example the :code:`mannings_n_at_sea_floor` can be set through the :code:`generate_linear_mannings_n()` method. Additionally, the user can set custom patches of values to the nodal attributes by passing a :code:`Polygon` or :code:`MultiPolygon` instance to the :code:`AdcircMesh.add_nodal_attributes_patch()` method.

``adcircpy.Fort15`` class
-------------------------
The :code:`Fort15` class meant to represent a ``fort.15`` file instance. As shown in the UML diagram in `Fort15 UML`_, this class contains attributes for all the configuration options for an  ADCIRC model run. Most of these options can be directly overridden by the users, with the exception of a very few. However, most parameters that cannot be overridden directly by the users provide methods that will ultimately change their value to what the user requires. One example of an option that is not directly configurable through the class attributes are the :code:`A00`, :code:`B00` and :code:`C00` factors. While the user cannot override any of these factors individually, the :code:`Fort15` class provides the :code:`set_time_weighting_factors_in_gcwe(A00, B00, C00)` which will make sure that the input parameters are consistent with what ADCIRC expects.

Another example of an important parameter that should not be directly overridden, but for which public methods are provided that will modify the options is the :code:`IM` parameter. This parameter is essentially a six digit integer that encodes multiple distinct configuration options. Additionally, this parameter has both single and double digit "shortcuts" to encode some of the most used options. These single and double digit codes can be expressed as six digit codes, therefore ADCIRCpy will only output the six digit version. This parameter is a very important one, and it might be difficult to memorize which digit does what, and `which digit position relates which option <https://wiki.adcirc.org/wiki/IM>`_. To deal with this, the :code:`Fort15` class attributes modify the final value of the :code:`IM` parameter, as configured by the users. To illustrate this, consider the case where the user instantiates the :code:`Fort15` class like this: :code:`fort15 = Fort15()`. Then the user can set the GWCE solver by setting :code:`fort15.gwce_solution_scheme = `implicit'`. By doing this setting, several parameters will change: the values of :code:`A00`, :code:`B00` and :code:`C00`, and the values of :code:`IM`. This particular parameter shows another example of how multiple lines in the ``fort.15`` need to change in order to keep consistency of the inputs. Other class attributes that affect the :code:`IM` parameter are: :code:`vertical_mode`, :code:`lateral_stress_in_gwce`, :code:`advection_in_gwce`, :code:`lateral_stress_in_momentum` and :code:`area_integration_in_momentum`, among others.

Although not recommended, if the user still wishes to bypass any sanity checks, all the ``fort.15`` parameters can be set directly by prefixing a single underscore in front of their corresponding Fort15 attribute. For example, the :code:`IM` parameter can be set directly by simply :code:`fort15._IM = 11111`. Note that using the single underscore prefix does not perform any sanity checks on the inputs, therefore by setting the values this way, ``fort.15`` output consistency, nor Fort15 functionality can be guaranteed, therefore it still advisable to use the provided public methods to set the configuration options for the ``fort.15`` file.

.. table:: GWCE solution scheme selection options

    ====================  ==== ==== ==== ======
    gwce.solution.scheme  A00  B00  C00   IM
    --------------------  ---- ---- ---- ------
    semi-implicit-legacy  0.35 0.3  0.35 511111
    semi-explicit         0    1    0    511112
    semi-implicit         0.5  0.5  0    511113
    ====================  ==== ==== ==== ======

.. _Fort15 UML:

.. figure:: figures/classes_adcircpy.fort15.png
  :alt: UML component diagram for the :code:`Fort15` class

  UML component diagram for the :code:`Fort15` class

``adcircpy.forcing`` module
---------------------------
The :code:`adcircpy.forcing` module contains a collection of modules, methods, functions, classes and procedures that directly relate to the different types of forcing that can be applied to the ADCIRC model. At the current state, the :code:`adcircpy.forcing` module contains two main modules that are of major interest to the users: the :code:`adcircpy.forcing.tides` module, and the :code:`adcircpy.forcing.winds` module. These will be discussed in sections `ADCIRCpy tides`_ and `ADCIRCpy winds`_ respectively.

.. _ADCIRCpy tides:

``adcircpy.forcing.tides`` module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :code:`adcircpy.forcing.tides` module contains one main public interface class :code:`adcircpy.forcing.tides.Tides` and two auxiliary classes for aggregating additional required tidal data: :code:`adcircpy.forcing.tides.HAMTIDE` (`ADCIRCpy HAMTIDE`_) and :code:`adcircpy.forcing.tides.TPXO` (`ADCIRCpy TPXO`_).

``adcircpy.forcing.tides.Tides`` class
""""""""""""""""""""""""""""""""""""""
The :code:`adcircpy.forcing.tides.Tides` class (abbreviated here as :code:`Tides`) provides the main interface for instantiating tidal elevation amplitudes and phases for a model run. This public class has methods such as :code:`use_constituent()`, :code:`use_all()` and :code:`use_major()` which can be used to enable the desired tidal constituents as forcings for the model. The :code:`Tides` class implements the computation of the initial conditions for tidal boundary. This class implements the same functionality as the ADCIRC ` :code:`tide_fac.f90` <https://www.dropbox.com/s/t2c1a6zo6tracs8/tide_fac.f>`_ program, which implements the equations of :cite:`Schureman1958` and it can be considered a direct "pythonic" port of this program. These equations provide the time dependent initial condition of the harmonic constituents of the tidal signal.

The :code:`Tides` class also implements the usage of helper classes (derived from the :code:`TidalDatabase` abstract class) in order to interpolate the harmonic elevation amplitude and phase for each tidal constituent at each boundary nodes. The concrete classes derived from :code:`TidalDatabase` class implement a :code:`get_amplitude(constituent, vertices)` and :code:`get_phase(constituent, vertices)` method that returns the sea surface elevation amplitude and phase for the ADCIRC model. The user can choose between the HAMTIDE database :cite:`Taguchi2014` or the TPXO database :cite:`Egbert2002` by using the :code:`database` keyword argument during instantiation of the :code:`Tides` class.

.. _ADCIRCpy HAMTIDE:

``adcircpy.forcing.tides.HAMTIDE`` class
""""""""""""""""""""""""""""""""""""""""
The :code:`adcircpy.forcing.tides.HAMTIDE` class is used by the :code:`Tides` class in order to interpolate harmonic constants to the boundary nodes, which are required for tidal model initialization. HAMTIDE database to obtain the information from the outputs of the Hamburg direct data Assimilation Methods for TIDEs model, which implements a methodology based on dynamical residuals as a component of the assimilation method to represent missing physics and their role in tidal dissipation by combination with measurements. This database includes amplitude and phase for sea-surface elevation and transport for the eight primary tidal constituents (N2, S2, M2, K2, P1, O1, K1, Q1). The data is provided as NetCDF files and through a live OpenDaP server.

.. _ADCIRCpy TPXO:

``adcircpy.forcing.tides.TPXO`` class
"""""""""""""""""""""""""""""""""""""
The :code:`adcircpy.forcing.tides.TPXO` provides amplitude and phase for sea-surface elevation and transport for the eight primary tidal constituents (N2, S2, M2, K2, P1, O1, K1, Q1), two long-period (Mf, Mm), and three non-linear (M4, MS4, MN4) harmonic components. The TPXO model uses a methodology based on least-squares best-fits of the Laplace Tidal Equation and with respect to observed altimetry data. The TPXO is a proprietary database that is provided in NetCDF and binary formats. The users are responsible for registering and provisioning a copy of the database for ADCIRCpy.

.. _ADCIRCpy winds:

``adcircpy.forcing.winds`` module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :code:`adcircpy.forcing.winds` module provides access to wind-related forcings for the model. In particular, the users might be interested in running a parametric wind field as a first approximation for the model and for stability checks. The class used to generate parametric wind fields in ADCIRC is described in `ADCIRCpy Best Track Forcing`_.

.. _ADCIRCpy Best Track Forcing:

``adcircpy.forcing.winds.BestTrackForcing`` class
"""""""""""""""""""""""""""""""""""""""""""""""""
The term "best track winds" refers to the set of spatially and temporally varying pressure and speed wind fields that are used to approximate planetary boundary layer vortices (hurricanes/typhoons). These best track winds are computed from sparse observational data. ADCIRCpy supports the generation of the input files required by ADCIRC in order to generate model runs using parametric wind models. The input files generated by ADCIRCpy uses by default the asymmetric wind field equations by setting the NWS parameter to the value of :code:`20` in the ``fort.15`` file.

``adcircpy.driver`` module
---------------------------
The :code:`adcircpy.driver` module consists of all methods and functions relating to ADCIRC forcing to  generate an overlying task submission framework and configuration directories. :code:`adcircpy.driver` contains module :code:`adcircpy.driver.AdcircRun`.

``adcircpy.driver.AdcircRun`` class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:code:`adcircpy.driver.AdcircRun` can determine the tidal constituents (phases and amplitudes) at either selected points in the domain for velocity and elevation or at all points in the domain for velocity and elevation.

``adcircpy.outputs`` module
---------------------------
:code:`adcircpy.outputs` module contain following output files :code:`OutputFactory`, :code:`OutputStations`, :code:`StationTimeseries`, :code:`ElevationStationsTimeseries`, :code:`HarmonicConstituentsElevationStations`, :code:`HarmonicConstituentsSurface`, :code:`Maxele`, :code:`ScalarSurfaceExtrema`, :code:`ScalarSurfaceTimeseries` and :code:`VelocityStations`. These files provide valuable information about the currents and tides at all points in a model domain. Regarding the output, ADCIRC gives excellent flexibility. For example, the user can define that water surface elevations and velocities be inscribed only at chosen points in the domain ``fort.61``. Alternatively, output at all positions in the environment can be requested ``fort.63`` for velocity and elevation output.

bibliography
------------

.. bibliography:: references.bib
