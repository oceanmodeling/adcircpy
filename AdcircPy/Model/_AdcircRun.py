import os
from datetime import datetime
import numpy as np
from scipy.interpolate import RectBivariateSpline
from netCDF4 import Dataset
from AdcircPy.Tides import TidalForcing

class _AdcircRun(object):
  def __init__(self, AdcircMesh, Tides=None, Winds=None, Waves=None, ElevationStationsOutput=None, VelocityStationsOutput=None, ElevationGlobalOutput=None, VelocityGlobalOutput=None, **kwargs):
    self.AdcircMesh = AdcircMesh
    self.init_Tides(Tides)
    self.init_Winds(Winds)
    # self.init_Waves(Waves)
    self.ElevationStationsOutput = ElevationStationsOutput
    self.init_fort15(**kwargs)

  def init_Tides(self, Tides):
    self.Tides = Tides
    self.init_TPXO()

  def init_fort15(self, **kwargs):
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
      self.IHOT=0
      self.write_fort15()

    with open(self.directory+'/fort.15.hotstart', 'w') as self.f:
      self.IHOT=567
      self.write_fort15()

  def init_TPXO(self):
    if self.AdcircMesh.ocean_boundaries is None or self.Tides is None:
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
          _x = (self.AdcircMesh.x[ocean_boundary]) % 360.
          _y = (self.AdcircMesh.y[ocean_boundary]) % 360.
          constituents[constituent]['ha'] = ha_interpolator.ev(_x, _y)
          constituents[constituent]['hp'] = hp_interpolator.ev(_x, _y)
      self.TPXO.append(constituents)

  def write_fort15(self):
    self.f.write('{:<32}! 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION\n'.format(self.RUNDES))
    self.f.write('{:<24}{:8}! 24 CHARACTER ALPANUMERIC RUN IDENTIFICATION\n'.format(self.RUNID,''))
    self.f.write('{:<32d}! NFOVER - NONFATAL ERROR OVERRIDE OPTION\n'.format(self.NFOVER))
    self.f.write('{:<32d}! NABOUT - ABREVIATED OUTPUT OPTION PARAMETER\n'.format(self.NABOUT))
    self.f.write('{:<32d}! NSCREEN - UNIT 6 OUTPUT OPTION PARAMETER\n'.format(self.NSCREEN))
    self.write_IHOT()
    self.f.write('{:<32d}! ICS - COORDINATE SYSTEM SELECTION PARAMETER\n'.format(self.ICS))
    self.f.write('{:<32d}! IM - MODEL SELECTION PARAMETER\n'.format(self.IM))
    self.f.write('{:<32d}! NOLIBF - BOTTOM FRICTION TERM SELECTION PARAM; before NWP==1, \'2\' was used\n'.format(self.NOLIBF))
    self.f.write('{:<32d}! NOLIFA - FINITE AMPLITUDE TERM SELECTION PARAMETER\n'.format(self.NOLIFA))
    self.f.write('{:<32d}! NOLICA - SPATIAL DERIVATIVE CONVECTIVE SELECTION PARAMETER\n'.format(self.NOLICA))
    self.f.write('{:<32d}! NOLICAT - TIME DERIVATIVE CONVECTIVE TERM SELECTION PARAMETER\n'.format(self.NOLICAT))
    self.write_NWP() # depends on fort.13 and a better implementation could be designed.
    self.f.write('{:<32d}! NCOR - VARIABLE CORIOLIS IN SPACE OPTION PARAMETER\n'.format(self.NCOR))
    self.write_NTIP()
    self.f.write('{:<32d}! NTIP - TIDAL POTENTIAL OPTION PARAMETER\n'.format(self.NTIP))
    self.write_NWS()
    self.write_NRAMP()
    self.f.write('{:<32.2f}! G - ACCELERATION DUE TO GRAVITY - DETERMINES UNITS\n'.format(self.G))
    self.write_TAU0()
    self.write_DTDP()
    self.f.write('{:<32.2f}! STATIM - STARTING TIME (IN DAYS)\n'.format(0))
    self.f.write('{:<32.2f}! REFTIM - REFERENCE TIME (IN DAYS)\n'.format(0))
    self.write_RNDAY()
    self.write_DRAMP()
    self.f.write('{:<4.2f} {:<4.2f} {:<4.2f} {:<17}! TIME WEIGHTING FACTORS FOR THE GWCE EQUATION\n'.format(self.A00, self.B00, self.C00, ''))
    self.write_H0_VELMIN()
    self.write_SLAM0_SFEA0()
    self.write_FFACTOR()
    self.write_ESLM()
    self.write_CORI()
    self.write_NTIF()
    self.write_NBFR()
    self.f.write('{:<32.1f}! ANGINN : INNER ANGLE THRESHOLD\n'.format(self.ANGINN))
    self.write_station_outputs()
    self.write_global_outputs()
    self.write_harmonic_outputs()
    self.write_hotstart_parameters()
    self.write_iteration_parameters()
    self.write_netcdf_parameters()
    self.write_fortran_namelists()
    self.f.write('\n')

  def write_IHOT(self):
    if self.IHOT==0:
      self.f.write('{:<32d}'.format(0))
    elif self.IHOT==567:
      self.f.write('{:<32d}'.format(567))
    self.f.write('! IHOT - HOT START PARAMETER\n')

  def write_NWP(self):
    if self.AdcircMesh.fort13 is None:
      self.f.write('{:<32d}'.format(0))
      self.f.write('! NWP - VARIABLE BOTTOM FRICTION AND LATERAL VISCOSITY OPTION PARAMETER; default 0\n')
    else:
      NWP = len(self.AdcircMesh.fort13.keys())
      self.f.write('{:<32d}'.format(NWP))
      self.f.write('! NWP - VARIABLE BOTTOM FRICTION AND LATERAL VISCOSITY OPTION PARAMETER; default 0\n')
      for attribute in self.AdcircMesh.fort13.keys():
        self.f.write('{}\n'.format(attribute))

  def write_NTIP(self):
    if self.NTIP is None:
      if self.Tides is None:
        self.NTIP=0
      else:
        self.NTIP=1
    elif self.NTIP=='fort.24':
      self.NTIP=2
    self.f.write('{:<32d}! NTIP - TIDAL POTENTIAL OPTION PARAMETER\n'.format(self.NTIP))

  def write_NWS(self):
    if self.IHOT==0 or self.Winds is None:
      self.f.write('{:<32d}'.format(0))
    else:
      self.f.write('{:<32}'.format('Please set manually.'))
    self.f.write('! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER\n')

  def write_NRAMP(self):
    if self.IHOT==0:
      if self.Tides is not None:
        self.NRAMP=1
       
    elif self.Winds is not None:
      self.NRAMP=8
      self.f.write('{:<32d}'.format(self.NRAMP))
    
    if self.NRAMP is not None:
      self.f.write('{:<32d}'.format(self.NRAMP))
    else:
      self.f.write('{:<32}'.format(''))
    self.f.write('! NRAMP - RAMP FUNCTION OPTION\n'.format(self.NRAMP))

  def write_TAU0(self):
    if self.AdcircMesh.fort13 is not None:
      if 'primitive_weighting_in_continuity_equation' in list(self.AdcircMesh.fort13.keys()):
        self.f.write('{:<32d}! TAU0 - WEIGHTING FACTOR IN GWCE; original, 0.005\n'.format(-3))
    elif self.TAU0==-5:
      self.f.write('{:<10.3f}'.format(self.Tau0FullDomainMin))
      self.f.write('{:<10.3f}\n'.format(self.Tau0FullDomainMax))
    else:
      self.f.write('{:<32}! TAU0 - WEIGHTING FACTOR IN GWCE; original, 0.005\n'.format(0.005))

  def write_DTDP(self):
    """
    The reason this is implemented separately is so that we can
    implement an optimal DTDP calculation in the future, based on
    the provided mesh.
    """
    self.f.write('{:<32.1f}! DT - TIME STEP (IN SECONDS)\n'.format(self.DTDP))

  def write_RNDAY(self):
    # Based on tides or based on winds? 
    # What if this is a met-only run without tides?
    if self.Tides is not None:
      if self.IHOT==0:
        RNDAY = (self.Tides.start_date - self.Tides.spinup_date).total_seconds()/(60*60*24)
      elif self.IHOT==567:
        RNDAY = (self.Tides.end_date - self.Tides.spinup_date).total_seconds()/(60*60*24)
      self.f.write('{:<32.2f}'.format(RNDAY))
    else:
      self.f.write('{:<32}'.format(''))
    self.f.write('! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)\n')

  def write_DRAMP(self):
    if self.Tides is not None:
      self.DRAMP = ((2/3)*(self.Tides.start_date - self.Tides.spinup_date).total_seconds())/(60*60*24)

    if self.NRAMP==1:
      self.f.write('{:<32.1f}'.format(self.DRAMP))
   
    elif self.NRAMP==8:
      if self.DUnRampMete is None:
        self.DUnRampMete = (self.Tides.start_date - self.Tides.spinup_date).days
      if self.DRAMPElev is None:
        self.DRAMPElev = self.DRAMP
      if self.DRAMPTip is None:
        self.DRAMPTip = self.DRAMP
      self.f.write('{:.1f} '.format(self.DRAMP))
      self.f.write('{:.1f} '.format(self.DRAMPExtFlux))
      self.f.write('{:.1f} '.format(self.FluxSettlingTime))
      self.f.write('{:.1f} '.format(self.DRAMPIntFlux))
      self.f.write('{:.1f} '.format(self.DRAMPElev))
      self.f.write('{:.1f} '.format(self.DRAMPTip))
      self.f.write('{:.1f} '.format(self.DRAMPMete))
      self.f.write('{:.1f} '.format(self.DRAMPWRad))
      self.f.write('{:.1f} '.format(self.DUnRampMete))
    else:
      self.f.write('{:<32}'.format(''))
    self.f.write('! DRAMP [DRAMPExtFlux, FluxSettlingTime,DRAMPIntFlux,DRAMPElev,DRAMPTip,DRAMPMete,DRAMPWRad,DRAMPUnMete] DURATION OF RAMP FUNCTION (IN DAYS)\n')

  def write_H0_VELMIN(self):
    if self.NOLIFA in [0,1]:
      self.f.write('{:<32.4f}'.format(self.H0))
    elif self.NOLIFA in [2,3]:
      self.f.write('{:<4.3f} 0 0 {:4.3f}{:<17}'.format(self.H0, self.VELMIN,''))
    self.f.write('! H0, NODEDRYMIN, NODEWETRMP, VELMIN\n')

  def write_SLAM0_SFEA0(self):
    # This is just the center of mass of the mesh.
    self.SLAM0 = np.mean(self.AdcircMesh.x)
    self.SFEA0 = np.mean(self.AdcircMesh.y)
    self.f.write('{:<4.1f} {:<4.1f}{:<22}'.format(self.SLAM0, self.SFEA0,''))
    self.f.write('! SLAM0,SFEA0 - CENTER OF CPP PROJECTION (NOT USED IF ICS=1, NTIP=0, NCOR=0)\n')

  def write_FFACTOR(self):
    if self.NOLIBF==0:
      self.f.write('{:<32.4f}'.format(self.TAU))
    elif self.NOLIBF==1:
      self.f.write('{:<32.4f}'.format(self.FFACTOR))
    elif self.NOLIBF==2:
      self.f.write('{:<6.4f} {:<6.4f} {:<6.4f} {:<6.4f} {:<3}'.format(self.FFACTOR, self.HBREAK, self.FTHETA, self.FGAMMA, ''))
    self.f.write('! FFACTOR\n')

  def write_ESLM(self):
    if self.IM in [0,1,2]:
      self.f.write('{:<32.2f}'.format(self.ESLM))
    elif self.IM in [10]:
      self.f.write('{:<2d} {}{:26}'.format(self.ESLM, self.ESLS,''))
    self.f.write('! ESL - LATERAL EDDY VISCOSITY COEFFICIENT; IGNORED IF NWP =1\n')

  def write_CORI(self):
    if self.NCOR==1:
      self.f.write('{:<32.1f}'.format(0))
    else:
      self.f.write('set parameter manually. ')
    self.f.write('! CORI - CORIOLIS PARAMETER - IGNORED IF NCOR = 1\n')

  def write_NTIF(self):
    """ 
    This segment is confusing and redundant.
    Could have added the extra parameters to NBFR or have them
    as part of the source code, since these are constants.
    """
    if self.NTIP==0:
      self.f.write('{:<32d}! NUMBER OF TIDAL POTENTIAL CONSTITUENTS BEING FORCED\n'.format(0))
    elif self.NTIP==1:
      NTIF=list()
      for constituent in self.Tides.keys():
        if 'earth_tidal_potential_reduction_factor' in self.Tides[constituent].keys():
          NTIF.append(constituent)
      self.f.write('{:<32d}! NUMBER OF TIDAL POTENTIAL CONSTITUENTS BEING FORCED\n'.format(len(NTIF)))
      for i, constituent in enumerate(NTIF):
        self.f.write('{:<32}'.format(constituent))
        if i == 0:
          self.f.write('! ALPHANUMERIC DESCRIPTION OF TIDAL POTENTIAL CONSTIT.\n')
        else:
          self.f.write('\n')
        self.f.write('{:>10.6f}'.format(self.Tides[constituent]['tidal_potential_amplitude']))
        self.f.write('{:>19.15f}'.format(self.Tides[constituent]['orbital_frequency']))
        self.f.write('{:>7.3f}'.format(self.Tides[constituent]['earth_tidal_potential_reduction_factor']))
        self.f.write('{:>9.5f}'.format(self.Tides[constituent]['nodal_factor']))
        self.f.write('{:>11.2f}'.format(self.Tides[constituent]['greenwich_term']))
        self.f.write('\n')
    elif NTIP==2 or NTIP=='fort.22':
      self.f.write('reading from fort.24, set parameter mannually.  ! NUMBER OF TIDAL POTENTIAL CONSTITUENTS BEING FORCED\n')

  def write_NBFR(self):
    if self.Tides is None:
      self.f.write('{:<32d}! NBFR: bnd forcing\n'.format(0))
    else:
      self.f.write('{:<32d}! NBFR: bnd forcing\n'.format(len(self.Tides.constituents)))
      for constituent in self.Tides.constituents:
        self.f.write('{}\n'.format(constituent))
        self.f.write('{:>19.15f}'.format(self.Tides[constituent]['orbital_frequency']))
        self.f.write('{:>9.5f}'.format(self.Tides[constituent]['nodal_factor']))
        self.f.write('{:>11.2f}'.format(self.Tides[constituent]['greenwich_term']))
        self.f.write('\n')
      for constituent in self.Tides.constituents:
        self.f.write('{}\n'.format(constituent.lower()))
        for boundary in self.TPXO:
          for i in range(boundary[constituent]['ha'].size):
            self.f.write('{:>11.6f}'.format(boundary[constituent]['ha'][i]))
            self.f.write('{:>14.6f}'.format(boundary[constituent]['hp'][i]))
            self.f.write('\n')

  def write_station_outputs(self):
    if self.ElevationStationsOutput is None:
      NOUTE=0
      TOUTSE=0
      TOUTFE=0
      NSPOOLE=0
      NSTAE=0
    elif self.ElevationStationsOutput is not None:
      if self.IHOT==0 or self.Tides is None:
        NOUTE=0
        TOUTSE=0
        TOUTFE=0
        NSPOOLE=0
        NSTAE=0
      else:
        if self.ElevationStationsOutput.netcdf==True:
          NOUTE=-5
        else:
          NOUTE=-1
        TOUTSE=(self.Tides.start_date - self.Tides.spinup_date).total_seconds()/(60*60*24)
        TOUTFE=(self.Tides.end_date - self.Tides.spinup_date).total_seconds()/(60*60*24)
        NSPOOLE=self.ElevationStationsOutput.sampling_frequency.seconds/self.DTDP
        NSTAE=len(self.ElevationStationsOutput.keys())
    self.f.write('{:<3d}'.format(NOUTE))
    self.f.write('{:<6.1f}'.format(TOUTSE))
    self.f.write('{:<8.2f}'.format(TOUTFE))
    self.f.write('{:<6.1f}'.format(NSPOOLE))
    self.f.write('{:<9}{}\n'.format('','! NOUTE,TOUTSE,TOUTFE,NSPOOLE:ELEV STATION OUTPUT INFO (UNIT 61)'))
    self.f.write('{:<32d}{}\n'.format(NSTAE, '! TOTAL NUMBER OF ELEVATION RECORDING STATIONS'))
    if self.ElevationStationsOutput is not None and self.IHOT==567:
      for station in self.ElevationStationsOutput.keys():
        self.f.write('{:<13.6f}'.format(self.ElevationStationsOutput[station]['x']))
        self.f.write('{:<13.6f}'.format(self.ElevationStationsOutput[station]['y']))
        self.f.write('! {}\n'.format(station))

  def write_global_outputs(self):
    pass

  def write_harmonic_outputs(self):
    pass

  def write_hotstart_parameters(self):
    pass

  def write_iteration_parameters(self):
    pass

  def write_netcdf_parameters(self):
    pass

  def write_fortran_namelists(self):
    pass
