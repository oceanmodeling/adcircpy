import os
import abc
from datetime import datetime, timedelta
from collections import Counter
import numpy as np
from scipy.interpolate import RectBivariateSpline, griddata, interp1d
from netCDF4 import Dataset
from AdcircPy.Model.ElevationStationsOutput import ElevationStationsOutput
from AdcircPy.Model.VelocityStationsOutput import VelocityStationsOutput
from AdcircPy.Model.MeteorologicalStationsOutput import MeteorologicalStationsOutput
from AdcircPy.Model.ElevationGlobalOutput import ElevationGlobalOutput
from AdcircPy.Model.VelocityGlobalOutput import VelocityGlobalOutput
from AdcircPy.Model.MeteorologicalGlobalOutput import MeteorologicalGlobalOutput
from AdcircPy.core import ServerConfiguration

class _AdcircRun(metaclass=abc.ABCMeta):
  def __init__(self, AdcircMesh, ElevationStationsOutput=None, VelocityStationsOutput=None, MeteorologicalStationsOutput=None,
                                 ElevationGlobalOutput=None, VelocityGlobalOutput=None, MeteorologicalGlobalOutput=None,
                                  netcdf=True, **kwargs):
    self.AdcircMesh = AdcircMesh
    self.ElevationStationsOutput = ElevationStationsOutput
    self.VelocityStationsOutput  = VelocityStationsOutput
    self.MeteorologicalStationsOutput  = MeteorologicalStationsOutput
    self.ElevationGlobalOutput   = ElevationGlobalOutput
    self.VelocityGlobalOutput    = VelocityGlobalOutput
    self.MeteorologicalGlobalOutput    = MeteorologicalGlobalOutput
    self.netcdf=netcdf
    self._init_fort15(**kwargs)
    self._init_NWS()
    self._init_DRAMP()
    self._init_TPXO()
  
  def dump(self, directory, printf=False):
    self.directory = directory
    os.makedirs(self.directory, exist_ok=True)
    
    with open(self.directory+'/fort.15.coldstart', 'w') as self.f:
      self.IHOT=0
      self._write_fort15()

    with open(self.directory+'/fort.15.hotstart', 'w') as self.f:
      self._set_IHOT()
      self._write_fort15()

    if hasattr(self, '_ServerConfiguration'):
      self._ServerConfiguration.dump(self.directory)
    
    if printf==True:
      with open(self.directory+'/fort.15.coldstart', 'r') as fort15:
        for line in fort15.read().splitlines():
          print(line)

      with open(self.directory+'/fort.15.hotstart', 'r') as fort15:
        for line in fort15.read().splitlines():
          print(line)

  def ServerConfiguration(self, *argv, **kwargs):
    """
    This method is used to provide server configuration parameters.
    More docstring will be provided in the future.
    """
    self._ServerConfiguration = ServerConfiguration(self, *argv, **kwargs)

  @abc.abstractmethod
  def _init_NWS(self):
    """
    self.NWS = ??
    """

  @abc.abstractmethod
  def _init_DRAMP(self):
    """
    self.NWS = ??
    """

  @abc.abstractmethod
  def _write_NRAMP(self):
    """
    Depends on whether wind forcing is present or not.
    """

  @abc.abstractmethod
  def _write_RNDAY(self):
    """
    Depends on forcings.
    """

  @abc.abstractmethod
  def _write_DRAMP(self):
    """
    Gets untangled on daughter classes.
    """
  def _init_fort15(self, **kwargs):
    self.RUNDES = kwargs.pop("RUNDES", self.AdcircMesh.description)  # UPTO 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION
    self.RUNID = kwargs.pop("RUNID", "generated on {}".format(datetime.now().strftime('%Y-%m-%d')))    # UPTO 24 CHARACTER ALPANUMERIC RUN IDENTIFICATION
    self.NFOVER = kwargs.pop("NFOVER", 1) # NFOVER - NONFATAL ERROR OVERRIDE OPTION
    self.WarnElev = kwargs.pop("WarnElev", None) # –DDEBUG_WARN_ELEV
    self.iWarnElevDump = kwargs.pop("iWarnElevDump", None)        # –DDEBUG_WARN_ELEV
    self.WarnElevDumpLimit = kwargs.pop("WarnElevDumpLimit", None) # –DDEBUG_WARN_ELEV
    self.ErrorElev = kwargs.pop("ErrorElev", None) # –DDEBUG_WARN_ELEV
    self.NABOUT = kwargs.pop("NABOUT", 1) # NABOUT - ABREVIATED OUTPUT OPTION PARAMETER
    self.NSCREEN = kwargs.pop("NSCREEN", 100) # NSCREEN - OUTPUT TO UNIT 6 PARAMETER
    self.IHOT = kwargs.pop("IHOT", None) # provided by package
    self.ICS = kwargs.pop("ICS", 2) # Coordinate system. 1 Cartesian, 2 spherical.
    self.IM = kwargs.pop("IM", 0) # Defaults to 2D Barotropic.
    self.IDEN = kwargs.pop("IDEN", None) # used only for 3D Baroclinic
    self.NOLIBF = kwargs.pop("NOLIBF", 1)  # bottom friction; 0:linear, 1:quadratic, 2:hybrid non-linear
    self.NOLIFA = kwargs.pop("NOLIFA", 2)  # wetting/drying; 0: wet/dry off no finite amplitude, 1: wet/dry off with finite amplitude, 2: wet/dry on. 
    self.NOLICA = kwargs.pop("NOLICA", 1)  # Advection; 0:OFF, 1:ON
    self.NOLICAT = kwargs.pop("NOLICAT", 1) # Advection time derivative; 0:OFF, 1:ON
    self.NWP = kwargs.pop("NWP", None)  # Not used, provided by package
    self.NCOR = kwargs.pop("NCOR", 1)  # Coriolis; 0:spatially constant, 1:spatially variable
    self.NTIP = kwargs.pop("NTIP", None)  # Tidal potentials and self attraction; 0:OFF, 1:TidalPotential, 2:use_fort.24
    self.NRAMP = kwargs.pop("NRAMP", None) # Provided by package. 0 for coldstart, 8 for hotstart
    self.G = kwargs.pop("G", 9.81) # G - ACCELERATION DUE TO GRAVITY - DETERMINES UNITS
    self.TAU0 = kwargs.pop("TAU0", None)  # Many options, need to set on package. It actually depends on a tau0_gen.f subroutine
    self.Tau0FullDomainMin = kwargs.pop("Tau0FullDomainMin", 0.005) # using suggested default
    self.Tau0FullDomainMax = kwargs.pop("Tau0FullDomainMax", 0.2) # using suggested default
    self.DTDP = kwargs.pop("DTDP", 2) # May be provided by package in the future.
    self.STATIM = kwargs.pop("STATIM", None) # Package provided on date initialization.
    self.REFTIM = kwargs.pop("REFTIM", None) # Package provided on date initialization.
    self.WTIMINC = kwargs.pop("WTIMINC", None) # Package provided during met forcing init.
    self.RSTIMINC = kwargs.pop("RSTIMINC", None) # Package provided during wave initialization
    self.RNDAY = kwargs.pop("RNDAY", None) # Provided by package on date initialization.
    self.DRAMP = kwargs.pop("DRAMP", None) # Provided by package 
    self.DRAMPExtFlux = kwargs.pop("DRAMPExtFlux", 0) # Ramp external boundary fluxes
    self.FluxSettlingTime = kwargs.pop("FluxSettlingTime", 0) # 
    self.DRAMPIntFlux = kwargs.pop("DRAMPIntFlux", 0) # Ramp internal fluxes
    self.DRAMPElev = kwargs.pop("DRAMPElev", None) # Ramps for boundary
    self.DRAMPTip = kwargs.pop("DRAMPTip", None) # Ramp for tidal potentials
    self.DRAMPMete = kwargs.pop("DRAMPMete", 1.) # ramp wind an pressure
    self.DRAMPWRad = kwargs.pop("DRAMPWRad", 0) # Ramp for waves
    self.DUnRampMete = kwargs.pop("DUnRampMete", None) # meteo ramp offset relative to coldstart
    self.A00 = kwargs.pop("A00", 0.35) # k+1 TIME WEIGHTING FACTORS FOR THE GWCE EQUATION
    self.B00 = kwargs.pop("B00", 0.30) # k   TIME WEIGHTING FACTORS FOR THE GWCE EQUATION
    self.C00 = kwargs.pop("C00", 0.35) # k-1 TIME WEIGHTING FACTORS FOR THE GWCE EQUATION
    self.H0 = kwargs.pop("H0", 0.05) # water level threshold for wet/dry
    self.VELMIN = kwargs.pop("VELMIN", 0.05) # velocity threshold for wetting, used for increasing model stability from wetting/dry
    self.SLAM0 = kwargs.pop("SLAM0", None) # longitude on which the CPP coordinate projection is centered (in degrees) if ICS = 2.
    self.SFEA0 = kwargs.pop("SFEA", None) # longitude on which the CPP coordinate projection is centered (in degrees) if ICS = 2.
    self.TAU = kwargs.pop("TAU", self.TAU0) # Not commonly used, setting recommended value.
    self.FFACTOR = kwargs.pop("FFACTOR", 0.0025) # replaces self.TAU in fort.15
    self.HBREAK = kwargs.pop("HBREAK", 1) # Not commonly used, setting recommended value.
    self.FTHETA = kwargs.pop("FTHETA", 10) # Not commonly used, setting recommended value.
    self.FGAMMA = kwargs.pop("FGAMMA", 1./3.) # Not commonly used, setting recommended value.
    self.ESLM = kwargs.pop("ESLM", 10.) # Horizontal eddy viscosity constant for momentum equation
    self.ESLS = kwargs.pop("ESLS", '???') # Horizontal eddy diffusivity constant for transport equation, used only if IM=10
    self.CORI = kwargs.pop("CORI", 0.0001)  # Coriolis constant for NCOR=0
    self.NTIF = kwargs.pop("NTIF", None) # Provided by module
    self.NBFR = kwargs.pop("NBFR", None) # Provided by module
    self.ANGINN = kwargs.pop("ANGINN", 110.) # Inner angle velocity threshold
    self.FMV = kwargs.pop("FMV", 0)
    self.NHAINC = kwargs.pop("NHAINC", None)
    self.NFFR = kwargs.pop("NFFR", None) # inflow boundaries forcing terms
    self.ITITER = kwargs.pop("ITITER", 1) # Solver type: -1:lumped, 1:ITPACKV2D
    self.ISLDIA = kwargs.pop("ISLDIA", 0) # fort.33 verbosity output level for solver (ITPACKV2D)
    self.CONVCR = kwargs.pop("CONVCR", .000001) # convergence criteria
    self.ITMAX = kwargs.pop("ITMAX", 25) # maximum number of iterations per timestep
    self.ILUMP = kwargs.pop("ILUMP", 0) # not documented.
    self.projtitle = kwargs.pop("projtitle", "projtitle")
    self.projinst = kwargs.pop("projinst", "projinst")
    self.projsrc = kwargs.pop("projsrc", "projsrc")
    self.projhist = kwargs.pop("projhist", "projhist")
    self.projref = kwargs.pop("projref", "projref")
    self.projcom = kwargs.pop("projcom", "projcom")
    self.projhost = kwargs.pop("projhost", "projhost")
    self.conv = kwargs.pop("conv", "conv")
    self.contac= kwargs.pop("contac", "contac")
    if len(kwargs.keys())>0:
      raise Exception('Provided unknown keywords: {}'.format(list(kwargs.keys())))

  def _init_TPXO(self):
    """
    Method to generate the TPXO interpolation on the open boundaries.
    """
    if self.AdcircMesh.ocean_boundaries is None:
      self.TPXO=None; return
    elif self.TidalForcing is None:
      self.TPXO=None; return
    else:
      nc = Dataset(self.TidalForcing.tpxo_path)
      tpxo_constituents = nc['con'][:].tostring().decode('UTF-8').split()
      x = nc['lon_z'][:].reshape(nc['lon_z'].size)
      y = nc['lat_z'][:].reshape(nc['lat_z'].size)
      self.TPXO=list()
      for ocean_boundary in self.AdcircMesh.ocean_boundaries:
        _x = (self.AdcircMesh.x[ocean_boundary]) % 360.
        _y = (self.AdcircMesh.y[ocean_boundary]) % 360.
        _idx = np.where(np.logical_and(
                         np.logical_and(np.min(_x)<=x, x<=np.max(_x)),
                         np.logical_and(np.min(_y)<=y, y<=np.max(_y))))
        constituents = dict()
        for constituent in self.TidalForcing.constituents:
          # There might be false mismatches between both databases. (tidal components who's names mismatch)      
          if constituent.lower() in tpxo_constituents:
            constituents[constituent] = dict()
            idx = tpxo_constituents.index(constituent.lower())
            ha = nc['ha'][idx,:,:].reshape(nc['ha'][idx,:,:].size)
            hp = nc['hp'][idx,:,:].reshape(nc['hp'][idx,:,:].size)
            constituents[constituent]['ha'] = griddata((x[_idx], y[_idx]), ha[_idx], (_x,_y), method='nearest')
            constituents[constituent]['hp'] = griddata((x[_idx], y[_idx]), hp[_idx], (_x,_y), method='nearest')
            # Rectilinear gives nicer results but tapers too strongly near shorelines and causes blowups.
            # ha_interpolator = RectBivariateSpline(x, y, nc['ha'][idx,:,:])
            # hp_interpolator = RectBivariateSpline(x, y, nc['hp'][idx,:,:])
            # constituents[constituent]['ha'] = ha_interpolator.ev(_x, _y)
            # constituents[constituent]['hp'] = hp_interpolator.ev(_x, _y)
        self.TPXO.append(constituents)

  def _write_fort15(self):
    self.f.write('{:<32}! 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION\n'.format(self.RUNDES))
    self.f.write('{:<24}{:8}! 24 CHARACTER ALPANUMERIC RUN IDENTIFICATION\n'.format(self.RUNID,''))
    self.f.write('{:<32d}! NFOVER - NONFATAL ERROR OVERRIDE OPTION\n'.format(self.NFOVER))
    self.f.write('{:<32d}! NABOUT - ABREVIATED OUTPUT OPTION PARAMETER\n'.format(self.NABOUT))
    self.f.write('{:<32d}! NSCREEN - UNIT 6 OUTPUT OPTION PARAMETER\n'.format(self.NSCREEN))
    self._write_IHOT()
    self.f.write('{:<32d}! ICS - COORDINATE SYSTEM SELECTION PARAMETER\n'.format(self.ICS))
    self.f.write('{:<32d}! IM - MODEL SELECTION PARAMETER\n'.format(self.IM))
    self.f.write('{:<32d}! NOLIBF - BOTTOM FRICTION TERM SELECTION PARAM; before NWP==1, \'2\' was used\n'.format(self.NOLIBF))
    self.f.write('{:<32d}! NOLIFA - FINITE AMPLITUDE TERM SELECTION PARAMETER\n'.format(self.NOLIFA))
    self.f.write('{:<32d}! NOLICA - SPATIAL DERIVATIVE CONVECTIVE SELECTION PARAMETER\n'.format(self.NOLICA))
    self.f.write('{:<32d}! NOLICAT - TIME DERIVATIVE CONVECTIVE TERM SELECTION PARAMETER\n'.format(self.NOLICAT))
    self._write_NWP() # depends on fort.13 and a better implementation could be designed.
    self.f.write('{:<32d}! NCOR - VARIABLE CORIOLIS IN SPACE OPTION PARAMETER\n'.format(self.NCOR))
    self._write_NTIP()
    self._write_NWS()
    self._write_NRAMP()
    self.f.write('! NRAMP - RAMP FUNCTION OPTION\n'.format(self.NRAMP))
    self.f.write('{:<32.2f}! G - ACCELERATION DUE TO GRAVITY - DETERMINES UNITS\n'.format(self.G))
    self._write_TAU0()
    self._write_DTDP()
    self.f.write('{:<32.2f}! STATIM - STARTING TIME (IN DAYS)\n'.format(0))
    self.f.write('{:<32.2f}! REFTIM - REFERENCE TIME (IN DAYS)\n'.format(0))
    self._write_WTIMINC()
    self._write_RNDAY()
    self.f.write('! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)\n')
    self._write_DRAMP()
    self.f.write('! DRAMP [DRAMPExtFlux, FluxSettlingTime,DRAMPIntFlux,DRAMPElev,DRAMPTip,DRAMPMete,DRAMPWRad,DRAMPUnMete] DURATION OF RAMP FUNCTION (IN DAYS)\n')
    self.f.write('{:<4.2f} {:<4.2f} {:<4.2f} {:<17}! TIME WEIGHTING FACTORS FOR THE GWCE EQUATION\n'.format(self.A00, self.B00, self.C00, ''))
    self._write_H0_VELMIN()
    self._write_SLAM0_SFEA0()
    self._write_FFACTOR()
    self._write_ESLM()
    self._write_CORI()
    self._write_NTIF()
    self._write_NBFR()
    self.f.write('{:<32.1f}! ANGINN : INNER ANGLE THRESHOLD\n'.format(self.ANGINN))
    self._write_ElevationStationsOutput()
    self._write_VelocityStationsOutput()
    self._write_ConcentrationStationsOutput()
    self._write_MeteorologicalStationsOutput()
    self._write_ElevationGlobalOutput()
    self._write_VelocityGlobalOutput()
    self._write_ConcentrationGlobalOutput()
    self._write_MeteorologicalGlobalOutput()
    self._write_HarmonicAnalysisOutputs()
    self._write_HotstartParams()
    self._write_iteration_parameters()
    self._write_netcdf_parameters()
    # self._write_fortran_namelists()
    self.f.write('\n')

  def _set_IHOT(self):
    if self.netcdf==True:
      self.IHOT=567
    else:
      self.IHOT=67
      
  def _write_IHOT(self):
    if self.IHOT==0:
      self.f.write('{:<32d}'.format(0))
    elif self.IHOT!=0:
      self.f.write('{:<32d}'.format(567))
    self.f.write('! IHOT - HOT START PARAMETER\n')

  def _write_NWP(self):
    if len(self.AdcircMesh.fort13.spinup_attributes) == 0 and \
        len(self.AdcircMesh.fort13.runtime_attributes) == 0:
      self.f.write('{:<32d}'.format(0))
      self.f.write('! NWP - VARIABLE BOTTOM FRICTION AND LATERAL VISCOSITY OPTION PARAMETER; default 0\n')
    else:
      if self.IHOT==0:
        NWP = len(self.AdcircMesh.fort13.spinup_attributes)
        self.f.write('{:<32d}'.format(NWP))
        self.f.write('! NWP - VARIABLE BOTTOM FRICTION AND LATERAL VISCOSITY OPTION PARAMETER; default 0\n')
        for attribute in self.AdcircMesh.fort13.spinup_attributes:
          self.f.write('{}\n'.format(attribute))
      else:
        NWP = len(self.AdcircMesh.fort13.runtime_attributes)
        self.f.write('{:<32d}'.format(NWP))
        self.f.write('! NWP - VARIABLE BOTTOM FRICTION AND LATERAL VISCOSITY OPTION PARAMETER; default 0\n')
        for attribute in self.AdcircMesh.fort13.runtime_attributes:
          self.f.write('{}\n'.format(attribute))

  def _write_NTIP(self):
    if self.NTIP is None:
      if self.TidalForcing is None:
        self.NTIP=0
      else:
        self.NTIP=1
    elif self.NTIP=='fort.24':
      self.NTIP=2
    self.f.write('{:<32d}! NTIP - TIDAL POTENTIAL OPTION PARAMETER\n'.format(self.NTIP))

  def _write_NWS(self):
    if self.IHOT==0:
      self.f.write('{:<32d}'.format(0))
    else:
      self.f.write('{:<32d}'.format(self.NWS))
    self.f.write('! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER\n')
  
  def _write_TAU0(self):
    def __check_primitive(attr):
      if 'primitive_weighting_in_continuity_equation' in attr:
        self.f.write('{:<32d}! TAU0 - WEIGHTING FACTOR IN GWCE; original, 0.005\n'.format(-3))
    if self.TAU0!=-5 and self.IHOT==0:
      if len(self.AdcircMesh.fort13.spinup_attributes)>0:
        __check_primitive(self.AdcircMesh.fort13.spinup_attributes)
    elif self.TAU0!=-5 and self.IHOT!=0:
      if len(self.AdcircMesh.fort13.runtime_attributes)>0:
        __check_primitive(self.AdcircMesh.fort13.runtime_attributes)
    elif self.TAU0==-5:
      self.f.write('{:<10.3f}'.format(self.Tau0FullDomainMin))
      self.f.write('{:<10.3f}\n'.format(self.Tau0FullDomainMax))
    else:
      self.f.write('{:<32}! TAU0 - WEIGHTING FACTOR IN GWCE; original, 0.005\n'.format(0.005))

  def _write_DTDP(self):
    """
    The reason this is implemented separately is so that we can
    implement an optimal DTDP calculation in the future, based on
    the provided mesh.
    """
    self.f.write('{:<32.1f}! DT - TIME STEP (IN SECONDS)\n'.format(self.DTDP))
  
  def _write_WTIMINC(self):
    if self.IHOT>0:
      if self.NWS==20:
        self.f.write('{:<5}'.format(self.fort22.start_date.year))
        self.f.write('{:<3}'.format(self.fort22.start_date.strftime('%m')))
        self.f.write('{:<3}'.format(self.fort22.start_date.strftime('%d')))
        self.f.write('{:<3}'.format(self.fort22.start_date.strftime('%H')))
        self.f.write('{:<2}'.format(1))
        self.f.write('{:<4}'.format(0.9))
        self.f.write('{:<1}'.format(1))
        self.f.write('{:<11}'.format(''))
      if self.NWS>0:
        self.f.write('! WTIMINC - meteorological data time increment, Geofactor=1 for NWS=20\n')

  def _write_H0_VELMIN(self):
    if self.NOLIFA in [0,1]:
      self.f.write('{:<32.4f}'.format(self.H0))
    elif self.NOLIFA in [2,3]:
      self.f.write('{:<4.3f} 0 0 {:4.3f}{:<17}'.format(self.H0, self.VELMIN,''))
    self.f.write('! H0, NODEDRYMIN, NODEWETRMP, VELMIN\n')

  def _write_SLAM0_SFEA0(self):
    # This is just the center of mass of the mesh.
    self.SLAM0 = np.mean(self.AdcircMesh.x)
    self.SFEA0 = np.mean(self.AdcircMesh.y)
    self.f.write('{:<4.1f} {:<4.1f}{:<22}'.format(self.SLAM0, self.SFEA0,''))
    self.f.write('! SLAM0,SFEA0 - CENTER OF CPP PROJECTION (NOT USED IF ICS=1, NTIP=0, NCOR=0)\n')

  def _write_FFACTOR(self):
    if self.NOLIBF==0:
      self.f.write('{:<32.4f}'.format(self.TAU))
    elif self.NOLIBF==1:
      self.f.write('{:<32.4f}'.format(self.FFACTOR))
    elif self.NOLIBF==2:
      self.f.write('{:<6.4f} {:<6.4f} {:<6.4f} {:<6.4f} {:<3}'.format(self.FFACTOR, self.HBREAK, self.FTHETA, self.FGAMMA, ''))
    self.f.write('! FFACTOR\n')

  def _write_ESLM(self):
    # Unclear whether ESLM should be float or integer.
    # All the examples show int but logically it looks like
    # it should be float... ??
    if self.IM in [0,1,2]:
      self.f.write('{:<32d}'.format(int(self.ESLM)))
    elif self.IM in [10]:
      self.f.write('{:<2d} {}{:26}'.format(self.ESLM, self.ESLS,''))
    self.f.write('! ESL - LATERAL EDDY VISCOSITY COEFFICIENT; IGNORED IF NWP =1\n')

  def _write_CORI(self):
    if self.NCOR==1:
      self.f.write('{:<32.1f}'.format(0))
    else:
      self.f.write('set parameter manually. ')
    self.f.write('! CORI - CORIOLIS PARAMETER - IGNORED IF NCOR = 1\n')

  def _write_NTIF(self):
    """ 
    This segment is confusing and redundant.
    Could have added the extra parameters to NBFR or have them
    as part of the source code, since these are constants.
    """
    if self.NTIP==0:
      self.f.write('{:<32d}! NUMBER OF TIDAL POTENTIAL CONSTITUENTS BEING FORCED\n'.format(0))
    elif self.NTIP==1:
      NTIF=list()
      for constituent in self.TidalForcing.keys():
        if 'earth_tidal_potential_reduction_factor' in self.TidalForcing[constituent].keys():
          NTIF.append(constituent)
      self.f.write('{:<32d}! NUMBER OF TIDAL POTENTIAL CONSTITUENTS BEING FORCED\n'.format(len(NTIF)))
      for i, constituent in enumerate(NTIF):
        self.f.write('{:<32}'.format(constituent))
        if i == 0:
          self.f.write('! ALPHANUMERIC DESCRIPTION OF TIDAL POTENTIAL CONSTIT.\n')
        else:
          self.f.write('\n')
        self.f.write('{:>10.6f}'.format(self.TidalForcing[constituent]['tidal_potential_amplitude']))
        self.f.write('{:>19.15f}'.format(self.TidalForcing[constituent]['orbital_frequency']))
        self.f.write('{:>7.3f}'.format(self.TidalForcing[constituent]['earth_tidal_potential_reduction_factor']))
        self.f.write('{:>9.5f}'.format(self.TidalForcing[constituent]['nodal_factor']))
        self.f.write('{:>11.2f}'.format(self.TidalForcing[constituent]['greenwich_term']))
        self.f.write('\n')
    elif NTIP==2 or NTIP=='fort.22':
      self.f.write('reading from fort.24, set parameter mannually.  ! NUMBER OF TIDAL POTENTIAL CONSTITUENTS BEING FORCED\n')

  def _write_NBFR(self):
    if self.TidalForcing is None:
      self.f.write('{:<32d}! NBFR: bnd forcing\n'.format(0))
    else:
      self.f.write('{:<32d}! NBFR: bnd forcing\n'.format(len(self.TidalForcing.constituents)))
      for constituent in self.TidalForcing.constituents:
        self.f.write('{}\n'.format(constituent))
        self.f.write('{:>19.15f}'.format(self.TidalForcing[constituent]['orbital_frequency']))
        self.f.write('{:>9.5f}'.format(self.TidalForcing[constituent]['nodal_factor']))
        self.f.write('{:>11.2f}'.format(self.TidalForcing[constituent]['greenwich_term']))
        self.f.write('\n')
      for constituent in self.TidalForcing.constituents:
        self.f.write('{}\n'.format(constituent.lower()))
        for boundary in self.TPXO:
          for i in range(boundary[constituent]['ha'].size):
            self.f.write('{:>11.6f}'.format(boundary[constituent]['ha'][i]))
            self.f.write('{:>14.6f}'.format(boundary[constituent]['hp'][i]))
            self.f.write('\n')
  
  def __remove_stations_not_in_domain(self):
    minlon, maxlon, minlat, maxlat = self.AdcircMesh.get_extent()
    for station in self.__StationsOutput.keys():
      x = self.__StationsOutput[station]['x']
      y = self.__StationsOutput[station]['y']
      if minlon<x<maxlon == False or minlat<y<maxlat==False:
        del self.__StationsOutput[station]   

  def __write_StationsOutput(self):
    # Should NSPOOL be float or int?
    # forcing int, but is there a case where it should be float?
    self.__remove_stations_not_in_domain()
    NSTA = len(self.__StationsOutput.keys())
    if  NSTA > 0:
      if self.TidalForcing is not None:
        if self.__StationsOutput.netcdf==True:
          NOUT=-5
        else:
          NOUT=-1
        if self.__StationsOutput.spinup==True and self.IHOT==0:
          TOUTS=0
          TOUTF=(self.TidalForcing.start_date-self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
          NSPOOL=self.__StationsOutput.sampling_frequency.seconds/self.DTDP
        elif self.__StationsOutput.spinup==False and self.IHOT==0:
          NOUT=0
          TOUTS=0
          TOUTF=0
          NSPOOL=0
          NSTA=0
        else:
          TOUTS=(self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
          TOUTF=(self.TidalForcing.end_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
          NSPOOL=self.__StationsOutput.sampling_frequency.seconds/self.DTDP 
    else:
      NOUT=0
      TOUTS=0
      TOUTF=0
      NSPOOL=0
    self.f.write('{:<3d}'.format(NOUT))
    self.f.write('{:<6.1f}'.format(TOUTS))
    self.f.write('{:<8.2f}'.format(TOUTF))
    self.f.write('{:<6d}'.format(int(NSPOOL)))
    self.f.write('{:<9}{}\n'.format('',self.__StationsOutput._comment1))
    self.f.write('{:<32d}{}\n'.format(NSTA, self.__StationsOutput._comment2))
    if NSTA>0:
      if self.IHOT!=0 or self.__StationsOutput.spinup==True:
        for station in self.__StationsOutput.keys():
          self.f.write('{:<13.6f}'.format(self.__StationsOutput[station]['x']))
          self.f.write('{:<13.6f}'.format(self.__StationsOutput[station]['y']))
          self.f.write('{:<6}! {}\n'.format('',station))
    del self.__StationsOutput
  
  def _write_ElevationStationsOutput(self):
    if self.ElevationStationsOutput is None:
      self.ElevationStationsOutput = ElevationStationsOutput()
    self.__StationsOutput=self.ElevationStationsOutput
    self.__write_StationsOutput()

  def _write_VelocityStationsOutput(self):
    if self.VelocityStationsOutput is None:
      self.VelocityStationsOutput = VelocityStationsOutput()
    self.__StationsOutput=self.VelocityStationsOutput
    self.__write_StationsOutput()

  def _write_ConcentrationStationsOutput(self):
    if self.IM==10:
      raise Exception('When IM=10 this line needs to be developed.')

  def _write_MeteorologicalStationsOutput(self):
    if self.IHOT>0:
      if self.NWS>0:
        if self.MeteorologicalStationsOutput is None:
          self.MeteorologicalStationsOutput = MeteorologicalStationsOutput()
        self.__StationsOutput=self.MeteorologicalStationsOutput
        self.__write_StationsOutput()

  def __write_GlobalOutputs(self):
    if self.TidalForcing is not None and self.__GlobalOutputs is not None:
      if self.__GlobalOutputs.netcdf==True:
        NOUTG=-5
      else:
        NOUTG=-1
      if self.__GlobalOutputs.spinup==True and self.IHOT==0:
        TOUTSG=0
        TOUTFG=(self.TidalForcing.start_date-self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
        NSPOOLG=self.__GlobalOutputs.sampling_frequency.seconds/self.DTDP
      elif self.__GlobalOutputs.spinup==False and self.IHOT==0:
        NOUTG=0
        TOUTSG=0
        TOUTFG=0
        NSPOOLG=0
      else:
        TOUTSG=(self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
        TOUTFG=(self.TidalForcing.end_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
        NSPOOLG=self.__GlobalOutputs.sampling_frequency.seconds/self.DTDP
    else:
      NOUTG=0
      TOUTSG=0
      TOUTFG=0
      NSPOOLG=0
    self.f.write('{:<3d}'.format(NOUTG))
    self.f.write('{:<6.1f}'.format(TOUTSG))
    self.f.write('{:<8.2f}'.format(TOUTFG))
    self.f.write('{:<6d}'.format(int(NSPOOLG)))
    del self.__GlobalOutputs

  def _write_ElevationGlobalOutput(self):
    self.__GlobalOutputs=self.ElevationGlobalOutput
    self.__write_GlobalOutputs()
    self.f.write('{:<9}{}'.format('','! NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE : GLOBAL ELEVATION OUTPUT INFO (UNIT  63)\n'))

  def _write_VelocityGlobalOutput(self):
    self.__GlobalOutputs = self.VelocityGlobalOutput
    self.__write_GlobalOutputs()
    self.f.write('{:<9}{}'.format('','! NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV : GLOBAL VELOCITY  OUTPUT INFO (UNIT  64)\n'))

  def _write_ConcentrationGlobalOutput(self):
    if self.IM==10:
      raise NotImplementedError('Concentration global outputs not yet implemented.')

  def _write_MeteorologicalGlobalOutput(self):
    if self.IHOT>0:
      if self.NWS>0:
        self.__GlobalOutputs = self.MeteorologicalGlobalOutput
        self.__write_GlobalOutputs()
        self.f.write('{:<9}{}'.format('','! NOUTGM,TOUTSGM,TOUTFGM,NSPOOLGM : GLOBAL METEOROLOGICAL  OUTPUT INFO (UNIT  64)\n'))
 
  def _write_HarmonicAnalysisOutputs(self):
    # Is NHAINC always an integer?
    if self.IHOT!=0 and self.TidalForcing is not None:

      if self.ElevationStationsOutput is not None:
        if self.ElevationStationsOutput.harmonic_analysis==True:
          if self.netcdf==True:
            NHASE = 5
          else:
            NHASE = 1
        else:
          NHASE = 0
      else:
        NHASE = 0

      if self.VelocityStationsOutput is not None:
        if self.VelocityStationsOutput.harmonic_analysis==True:
          if self.netcdf==True:
            NHASV = 5
          else:
            NHASV = 1
        else:
          NHASV = 0
      else:
        NHASV = 0

      if self.ElevationGlobalOutput is not None:
        if self.ElevationGlobalOutput.harmonic_analysis==True:
          if self.netcdf==True:
            NHAGE = 5
          else:
            NHAGE = 1
        else:
          NHAGE = 0
      else:
        NHAGE = 0

      if self.VelocityGlobalOutput is not None:
        if self.VelocityGlobalOutput.harmonic_analysis==True:
          if self.netcdf==True:
            NHAGV = 5
          else:
            NHAGV = 1
        else:
          NHAGV = 0
      else:
        NHAGV = 0
    else:
      NHASE=0
      NHASV=0
      NHAGE=0
      NHAGV=0

    if 1 in [NHASE, NHASV, NHAGE, NHAGV] or 5 in [NHASE, NHASV, NHAGE, NHAGV]:
      NFREQ=len(self.TidalForcing.keys())
      THAS=(self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds() / (60*60*24)
      THAF=(self.TidalForcing.end_date - self.TidalForcing.spinup_date).total_seconds() / (60*60*24)
      if self.NHAINC is None:
        NHAINC=timedelta(minutes=6).seconds/self.DTDP

    else:
      NFREQ=0
      THAS=0
      THAF=0
      if self.NHAINC is None:
        NHAINC=0
      else:
        NHAINC=self.NHAINC

    self.f.write('{:<32d}'.format(NFREQ))
    self.f.write('! NFREQ - NUMBER OF CONSTITUENTS TO BE INCLUDED IN THE HARMONIC ANALYSIS OUTPUTS\n')
    if NFREQ>0:
      for constituent in self.TidalForcing.constituents:
        self.f.write('{}\n'.format(constituent))
        self.f.write('{:>19.15f}'.format(self.TidalForcing[constituent]['orbital_frequency']))
        self.f.write('{:>9.5f}'.format(self.TidalForcing[constituent]['nodal_factor']))
        self.f.write('{:>11.2f}'.format(self.TidalForcing[constituent]['greenwich_term']))
        self.f.write('\n')
    self.f.write('{:<6.1f}'.format(THAS))
    self.f.write('{:<6.1f}'.format(THAF))
    self.f.write('{:<6d}'.format(int(NHAINC)))
    self.f.write('{:<3.1f}'.format(self.FMV))
    self.f.write('{:<11}{}'.format('','! THAS,THAF,NHAINC,FMV - HARMONIC ANALYSIS PARAMETERS\n'))
    self.f.write('{:<4d}'.format(NHASE))
    self.f.write('{:<4d}'.format(NHASV))
    self.f.write('{:<4d}'.format(NHAGE))
    self.f.write('{:<4d}'.format(NHAGV))
    self.f.write('{:16}{}'.format('','! NHASE,NHASV,NHAGE,NHAGV - CONTROL HARMONIC ANALYSIS AND OUTPUT TO UNITS 51,52,53,54\n'))

  def _write_HotstartParams(self):
    if self.IHOT==0:
      if self.netcdf==True:
        NHSTAR=5
      else:
        NHSTAR=1
      if self.TidalForcing is not None:
        NHSINC=(self.TidalForcing.start_date-self.TidalForcing.spinup_date).total_seconds()/self.DTDP
      else:
        NHSINC=0
    else:
      NHSTAR=0
      NHSINC=0
    self.f.write('{:<4d}'.format(NHSTAR))
    self.f.write('{:<10d}'.format(int(np.around(NHSINC,0))))
    self.f.write('{:<18}{}'.format('','! NHSTAR,NHSINC - HOT START FILE GENERATION PARAMETERS\n'))
    
  def _write_iteration_parameters(self):
    self.f.write('{:<4d}'.format(self.ITITER))
    self.f.write('{:<4d}'.format(self.ISLDIA))
    self.f.write('{:<7.0E}'.format(self.CONVCR))
    self.f.write('{:<4d}'.format(self.ITMAX))
    self.f.write('{:<4d}'.format(self.ILUMP))
    self.f.write('{:<9}{}'.format('','! ITITER, ISLDIA, CONVCR, ITMAX, ILUMP - ALGEBRAIC SOLUTION PARAMETERS\n'))

  def _write_netcdf_parameters(self):
    if self.netcdf==True:
      self.f.write('{}\n'.format(self.projtitle))
      self.f.write('{}\n'.format(self.projinst))
      self.f.write('{}\n'.format(self.projsrc))
      self.f.write('{}\n'.format(self.projhist))
      self.f.write('{}\n'.format(self.projref))
      self.f.write('{}\n'.format(self.projcom))
      self.f.write('{}\n'.format(self.projhost))
      self.f.write('{}\n'.format(self.conv))
      self.f.write('{}\n'.format(self.contac))
    if self.TidalForcing is not None:
      self.f.write('{:<23}\n'.format(self.TidalForcing.spinup_date.strftime('%Y-%m-%d %H:%M:%S UTC')))

  def _write_fortran_namelists(self):
    pass
