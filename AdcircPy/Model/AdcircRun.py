import os
from datetime import datetime
import numpy as np
from scipy.interpolate import RectBivariateSpline
from netCDF4 import Dataset
from AdcircPy.Tides import TidalForcing

def init_Tides(self, Tides):
  self.Tides = Tides
  self._init_TPXO()

def init_fort15(self, **kwargs):
  self.RUNDES = kwargs.pop("RUNDES", self.AdcircMesh.description)  # UPTO 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION
  self.RUNID = kwargs.pop("RUNID", "generated on {}".format(datetime.now().strftime('%Y-%m-%d')))    # UPTO 24 CHARACTER ALPANUMERIC RUN IDENTIFICATION
  self.NFOVER = kwargs.pop("NFOVER", 1) # NFOVER - NONFATAL ERROR OVERRIDE OPTION
  self.WarnElev = kwargs.pop("WarnElev", None) # –DDEBUG_WARN_ELEV
  self.iWarnElevDump = kwargs.pop("iWarnElevDump", None)        # –DDEBUG_WARN_ELEV
  self.WarnElevDumpLimit = kwargs.pop("WarnElevDumpLimit", None) # –DDEBUG_WARN_ELEV
  self.ErrorElev = kwargs.pop("ErrorElev", None) # –DDEBUG_WARN_ELEV
  self.NABOUT = kwargs.pop("NABOUT", 1) # NABOUT - ABREVIATED OUTPUT OPTION PARAMETER
  self.NSCREEN = kwargs.pop("NSCREEN", 1) # NSCREEN - OUTPUT TO UNIT 6 PARAMETER
  self.IHOT = kwargs.pop("IHOT", None) # provided by package
  self.ICS = kwargs.pop("ICS", 2) # Coordinate system. 1 Cartesian, 2 spherical.
  self.IM = kwargs.pop("IM", 0) # Defaults to 2D Barotropic.
  self.IDEN = kwargs.pop("IDEN", None) # used only for 3D Baroclinic
  self.NOLIBF = kwargs.pop("NOLIBF", 2)  # bottom friction; 0:linear, 1:quadratic, 2:hybrid non-linear
  self.NOLIFA = kwargs.pop("NOLIFA", 2)  # wetting/drying; 0: wet/dry off no finite amplitude, 1: wet/dry off with finite amplitude, 2: wet/dry on. 
  self.NOLICA = kwargs.pop("NOLICA", 1)  # Advection; 0:OFF, 1:ON
  self.NOLICAT = kwargs.pop("NOLICAT", 1) # Advection time derivative; 0:OFF, 1:ON
  self.NWP = kwargs.pop("NWP", None)  # Not used, provided by package
  self.NCOR = kwargs.pop("NCOR", 1)  # Coriolis; 0:spatially constant, 1:spatially variable
  self.NTIP = kwargs.pop("NTIP", 1)  # Tidal potentials and self attraction; 0:OFF, 1:TidalPotential, 2:use_fort.24
  self.NWS = kwargs.pop("NWS", None) # Meteo forcing selection. Provided by package.
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
  self.DRAMPExtFlux = kwargs.pop("DRAMPExtFlux", None) # Ramp external boundary fluxes
  self.FluxSettlingTime = kwargs.pop("FluxSettlingTime", None) # 
  self.DRAMPIntFlux = kwargs.pop("DRAMPIntFlux", None) # Ramp internal fluxes
  self.DRAMPElev = kwargs.pop("DRAMPElev", None) # Ramps for boundary
  self.DRAMPTip = kwargs.pop("DRAMPTip", None) # Ramp for tidal potentials
  self.DRAMPMete = kwargs.pop("DRAMPMete", None) # ramp wind an pressure
  self.DRAMPWRad = kwargs.pop("DRAMPWRad", None) # Ramp for waves
  self.DUnRampMete = kwargs.pop("DUnRampMete", None) # meteo ramp offset relative to coldstart
  self.A00 = kwargs.pop("A00", 0.35) # k+1 TIME WEIGHTING FACTORS FOR THE GWCE EQUATION
  self.B00 = kwargs.pop("B00", 0.30) # k   TIME WEIGHTING FACTORS FOR THE GWCE EQUATION
  self.C00 = kwargs.pop("C00", 0.35) # k-1 TIME WEIGHTING FACTORS FOR THE GWCE EQUATION
  self.H0 = kwargs.pop("H0", 0.05) # water level threshold for wet/dry
  self.VELMIN = kwargs.pop("VELMIN", 0.05) # velocity threshold for wetting, used for increasing model stability from wetting/dry
  self.SLAM0 = kwargs.pop("SLAM0", -71.0) # longitude on which the CPP coordinate projection is centered (in degrees) if ICS = 2.
  self.SFEA0 = kwargs.pop("SFEA", 27.0) # longitude on which the CPP coordinate projection is centered (in degrees) if ICS = 2.
  self.TAU = kwargs.pop("TAU", self.TAU0) # Not commonly used, setting recommended value.
  self.CF = kwargs.pop("CF", 0.0025) # replaces self.TAU in fort.15
  self.HBREAK = kwargs.pop("HBREAK", 1) # Not commonly used, setting recommended value.
  self.FTHETA = kwargs.pop("FTHETA", 10) # Not commonly used, setting recommended value.
  self.FGAMMA = kwargs.pop("FGAMMA", 1./3.) # Not commonly used, setting recommended value.
  self.ESLM = kwargs.pop("ESLM", 10.) # Horizontal eddy viscosity constant for momentum equation
  self.ESLM = kwargs.pop("ESLS", None) # Horizontal eddy diffusivity constant for transport equation, used only if IM=10
  self.CORI = kwargs.pop("CORI", 0.0001)  # Coriolis constant for NCOR=0
  self.NTIF = kwargs.pop("NTIF", None) # Provided by module
  self.NBFR = kwargs.pop("NBFR", None) # Provided by module
  self.ANGINN = kwargs.pop("ANGINN", 110.) # Inner angle velocity threshold
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
  self.NCDATE = kwargs.pop("NCDATE", None) # set by package on datetime init

def dump(self, directory):
  self.directory = directory
  os.makedirs(directory, exist_ok=True)
  
  with open(self.directory+'/fort.15.coldstart', 'w') as self.f:
    self.run_type='coldstart'
    self._write_fort15()

  with open(self.directory+'/fort.15.hotstart', 'w') as self.f:
    self.run_type='hotstart'
    self._write_fort15()


def _init_TPXO(self):
  if self.AdcircMesh.ocean_boundaries is None:
    self.TPXO=None; return
  nc = Dataset(self.Tides.tpxo_path)
  tpxo_constituents = nc['con'][:].tostring().decode('UTF-8').split()
  x = nc['lon_z'][:]
  y = nc['lat_z'][:]
  x = np.linspace(np.min(x),np.max(x),num=x.shape[0])
  y = np.linspace(np.min(y),np.max(y),num=y.shape[1])
  self.TPXO=list()
  for ocean_boundary in self.AdcircMesh.ocean_boundaries:
    constituents = dict()
    for constituent in self.Tides.constituents:
      # There might be false mismatches between both databases. (tidal components who's names mismatch)      
      if constituent.lower() in tpxo_constituents:
        constituents[constituent] = dict()
        idx = tpxo_constituents.index(constituent.lower())
        ha_interpolator = RectBivariateSpline(x, y, nc['ha'][idx,:,:])
        hp_interpolator = RectBivariateSpline(x, y, nc['hp'][idx,:,:])
        constituents[constituent]['ha'] = ha_interpolator.ev(self.AdcircMesh.x[ocean_boundary], self.AdcircMesh.y[ocean_boundary])
        constituents[constituent]['hp'] = hp_interpolator.ev(self.AdcircMesh.x[ocean_boundary], self.AdcircMesh.y[ocean_boundary])
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
  self.f.write('{:<32d}! NTIP - TIDAL POTENTIAL OPTION PARAMETER\n'.format(self.NTIP))
  self._write_NWS()
  self._write_NRAMP()
  self.f.write('{:<32.2f}! G - ACCELERATION DUE TO GRAVITY - DETERMINES UNITS\n'.format(self.G))
  self._write_TAU0()
  self._write_DTDP()
  self.f.write('{:<32.2f}! STATIM - STARTING TIME (IN DAYS)\n'.format(0))
  self.f.write('{:<32.2f}! REFTIM - REFERENCE TIME (IN DAYS)\n'.format(0))
  self._write_RNDAY()
  self.f.write('\n')
  

def _write_IHOT(self):
  if self.run_type=='coldstart':
    self.f.write('{:<32d}'.format(0))
  elif self.run_type=='hotstart':
    self.f.write('{:<32d}'.format(567))
  self.f.write('! IHOT - HOT START PARAMETER\n')

def _write_NWP(self):
  if self.AdcircMesh.fort13 is None:
    self.f.write('{:<32d}'.format(0))
    self.f.write('! NWP - VARIABLE BOTTOM FRICTION AND LATERAL VISCOSITY OPTION PARAMETER; default 0\n')
  else:
    NWP = len(self.AdcircMesh.fort13.keys())
    self.f.write('{:<32d}'.format(NWP))
    self.f.write('! NWP - VARIABLE BOTTOM FRICTION AND LATERAL VISCOSITY OPTION PARAMETER; default 0\n')
    for attribute in self.AdcircMesh.fort13.keys():
      self.f.write('{}\n'.format(attribute))

def _write_NWS(self):
  if self.run_type=='coldstart' or self.Winds is None:
    self.f.write('{:<32d}'.format(0))
  else:
    self.f.write('{:<32}'.format('Please set manually.'))
  self.f.write('! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER\n')

def _write_NRAMP(self):
  if self.run_type=='coldstart' or self.Winds is None:
    self.f.write('{:<32d}'.format(1))
  else:
    self.f.write('{:<32d}'.format(8))
  self.f.write('! NRAMP - RAMP FUNCTION OPTION\n')

def _write_TAU0(self):
  if 'primitive_weighting_in_continuity_equation' in list(self.AdcircMesh.fort13.keys()):
    self.f.write('{:<32d}! TAU0 - WEIGHTING FACTOR IN GWCE; original, 0.005\n'.format(-3))
  elif self.TAU0==-5:
    self.f.write('{:<10.3f}'.format(self.Tau0FullDomainMin))
    self.f.write('{:<10.3f}\n'.format(self.Tau0FullDomainMax))
  else:
    self.f.write('{:<32}! TAU0 - WEIGHTING FACTOR IN GWCE; original, 0.005'.format('Please set manually.'))

def _write_DTDP(self):
  """
  The reason this is implemented separately is so that we can
  implement an optimal DTDP calculation in the future, based on
  the provided mesh.
  """
  self.f.write('{:<32.1f}! DT - TIME STEP (IN SECONDS)\n'.format(self.DTDP))

def _write_RNDAY(self):
  if self.Tides is not None:
    if self.run_type=="coldstart":
      RNDAY = (self.Tides.start_date - self.Tides.spinup_date).days
    elif self.run_type=="hotstart":
      RNDAY = (self.Tides.end_date - self.Tides.spinup_date).days
    self.f.write('{:<32.2f}'.format(RNDAY))
  
  elif self.Winds is not None:
    # May be a Met-Only run
    self.f.write('{:<32}'.format('Set manually for met-only run'))
  self.f.write('! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)\n')

  