# global imports
import numpy as np
from datetime import datetime

# local imports
from AdcircPy import Model

# unittest import
import os
from datetime import timedelta
import unittest


class _Fort15(object):
    """
    Note from the author:
        There are details about how ADCIRC parses the fort.15 file that were
        unknown to me at the time of writting. There might be formatting styles
        that need to be modified. Scientific notation was used wherever
        practical in order to maintain the largest accuracy possible. The
        decimal places for the scientific notations may require further
        optimization.
    """

    def __init__(self,
                 RUNDES=None,
                 RUNID=None,
                 NFOVER=0,
                 WarnElev=None,
                 iWarnElevDump=None,
                 WarnElevDumpLimit=None,
                 ErrorElev=None,
                 NABOUT=1,
                 NSCREEN=500,
                 ICS=None,
                 IM=0,
                 IDEN=None,
                 NOLIBF=None,
                 NOLIFA=2,
                 NOLICA=1,
                 NOLICAT=1,
                 NCOR=None,
                 NTIP=1,
                 G=9.81,
                 DTDP=2.,
                 STATIM=None,
                 REFTIM=None,
                 DRAMP=None,
                 DRAMPExtFlux=0.,
                 FluxSettlingTime=0.,
                 DRAMPIntFlux=0.,
                 DRAMPElev=0.,
                 DRAMPTip=0.,
                 DRAMPMete=0.,
                 DRAMPWRad=0.,
                 DUnRampMete=None,
                 A00=0.35,
                 B00=0.30,
                 C00=0.35,
                 H0=0.05,
                 NODEDRYMIN=10,
                 NODEWETRMP=10,
                 VELMIN=0.05,
                 SLAM0=None,
                 SFEA0=None,
                 FFACTOR=0.,
                 HBREAK=1,
                 FTHETA=10,
                 FGAMMA=1./3.,
                 ESLM=10.,
                 ESLC=None,
                 CORI=None,
                 NTIF='all',
                 ANGINN=110.,
                 THAS=None,
                 THAF=None,
                 NHAINC=None,
                 FMV=0,
                 NHSTAR=None,
                 NHSINC=None,
                 ITITER=1,
                 ISLDIA=0,
                 CONVCR=1.e-6,
                 ITMAX=25,
                 ILUMP=0,
                 NCPROJ=None,
                 NCINST=None,
                 NCSOUR=None,
                 NCHIST=None,
                 NCREF=None,
                 NCCOM=None,
                 NCHOST=None,
                 NCCONV=None,
                 NCCONT=None,
                 NCDATE=None,
                 FortranNamelists=None):
        self.RUNDES = RUNDES
        self.RUNID = RUNID
        self.NFOVER = NFOVER
        self.WarnElev = WarnElev
        self.iWarnElevDump = iWarnElevDump
        self.WarnElevDumpLimit = WarnElevDumpLimit
        self.ErrorElev = ErrorElev
        self.NABOUT = NABOUT
        self.NSCREEN = NSCREEN
        self.ICS = ICS
        self.IM = IM
        self.IDEN = IDEN
        self.NOLIBF = NOLIBF
        self.NOLIFA = NOLIFA
        self.NOLICA = NOLICA
        self.NOLICAT = NOLICAT
        self.NCOR = NCOR
        self.NTIP = NTIP
        self.G = G
        self.DTDP = DTDP
        self.STATIM = STATIM
        self.REFTIM = REFTIM
        self.DRAMP = DRAMP
        self.DRAMPExtFlux = DRAMPExtFlux
        self.FluxSettlingTime = FluxSettlingTime
        self.DRAMPIntFlux = DRAMPIntFlux
        self.DRAMPElev = DRAMPElev
        self.DRAMPTip = DRAMPTip
        self.DRAMPMete = DRAMPMete
        self.DRAMPWRad = DRAMPWRad
        self.DUnRampMete = DUnRampMete
        self.A00 = A00
        self.B00 = B00
        self.C00 = C00
        self.H0 = H0
        self.NODEDRYMIN = NODEDRYMIN
        self.NODEWETRMP = NODEWETRMP
        self.VELMIN = VELMIN
        self.SLAM0 = SLAM0
        self.SFEA0 = SFEA0
        self.FFACTOR = FFACTOR
        self.HBREAK = HBREAK
        self.FTHETA = FTHETA
        self.FGAMMA = FGAMMA
        self.ESLM = ESLM
        self.ESLC = ESLC
        self.CORI = CORI
        self.NTIF = NTIF
        self.ANGINN = ANGINN
        self.THAS = THAS
        self.THAF = THAF
        self.NHAINC = NHAINC
        self.FMV = FMV
        self.NHSTAR = NHSTAR
        self.NHSINC = NHSINC
        self.ITITER = ITITER
        self.ISLDIA = ISLDIA
        self.CONVCR = CONVCR
        self.ITMAX = ITMAX
        self.ILUMP = ILUMP
        self.NCPROJ = NCPROJ
        self.NCINST = NCINST
        self.NCSOUR = NCSOUR
        self.NCHIST = NCHIST
        self.NCREF = NCREF
        self.NCCOM = NCCOM
        self.NCHOST = NCHOST
        self.NCCONV = NCCONV
        self.NCCONT = NCCONT
        self.NCDATE = NCDATE

    @property
    def RUNDES(self):
        return self.__RUNDES

    @property
    def RUNID(self):
        return self.__RUNID

    @property
    def NFOVER(self):
        return self.__NFOVER

    @property
    def WarnElev(self):
        return self.__WarnElev

    @property
    def iWarnElevDump(self):
        return self.__iWarnElevDump

    @property
    def WarnElevDumpLimit(self):
        return self.__WarnElevDumpLimit

    @property
    def ErrorElev(self):
        return self.__ErrorElev

    @property
    def NABOUT(self):
        return self.__NABOUT

    @property
    def NSCREEN(self):
        return self.__NSCREEN

    @property
    def ICS(self):
        return self.__ICS

    @property
    def IM(self):
        return self.__IM

    @property
    def IDEN(self):
        return self.__IDEN

    @property
    def NOLIBF(self):
        return self.__NOLIBF

    @property
    def NOLIFA(self):
        return self.__NOLIFA

    @property
    def NOLICA(self):
        return self.__NOLICA

    @property
    def NOLICAT(self):
        return self.__NOLICAT

    @property
    def NCOR(self):
        return self.__NCOR

    @property
    def NTIP(self):
        return self.__NTIP

    @property
    def G(self):
        return self.__G

    @property
    def DTDP(self):
        return self.__DTDP

    @property
    def STATIM(self):
        return self.__STATIM

    @property
    def REFTIM(self):
        return self.__REFTIM

    @property
    def DRAMP(self):
        return self.__DRAMP

    @property
    def DRAMPExtFlux(self):
        return self.__DRAMPExtFlux

    @property
    def FluxSettlingTime(self):
        return self.__FluxSettlingTime

    @property
    def DRAMPIntFlux(self):
        return self.__DRAMPIntFlux

    @property
    def DRAMPElev(self):
        return self.__DRAMPElev

    @property
    def DRAMPTip(self):
        return self.__DRAMPTip

    @property
    def DRAMPMete(self):
        return self.__DRAMPMete

    @property
    def DRAMPWRad(self):
        return self.__DRAMPWRad

    @property
    def DUnRampMete(self):
        return self.__DUnRampMete

    @property
    def A00(self):
        return self.__A00

    @property
    def B00(self):
        return self.__B00

    @property
    def C00(self):
        return self.__C00

    @property
    def H0(self):
        return self.__H0

    @property
    def NODEDRYMIN(self):
        return self.__NODEDRYMIN

    @property
    def NODEWETRMP(self):
        return self.__NODEWETRMP

    @property
    def VELMIN(self):
        return self.__VELMIN

    @property
    def SLAM0(self):
        return self.__SLAM0

    @property
    def SFEA0(self):
        return self.__SFEA0

    @property
    def FFACTOR(self):
        return self.__FFACTOR

    @property
    def TAU(self):
        return self.__TAU

    @property
    def HBREAK(self):
        return self.__HBREAK

    @property
    def FTHETA(self):
        return self.__FTHETA

    @property
    def FGAMMA(self):
        return self.__FGAMMA

    @property
    def ESLM(self):
        return self.__ESLM

    @property
    def ESLC(self):
        return self.__ESLC

    @property
    def CORI(self):
        return self.__CORI

    @property
    def NTIF(self):
        return self.__NTIF

    @property
    def ANGINN(self):
        return self.__ANGINN

    @property
    def THAS(self):
        return self.__THAS

    @property
    def THAF(self):
        return self.__THAF

    @property
    def NHAINC(self):
        return self.__NHAINC

    @property
    def FMV(self):
        return self.__FMV

    @property
    def NHSTAR(self):
        return self.__NHSTAR

    @property
    def NHSINC(self):
        return self.__NHSINC

    @property
    def ITITER(self):
        return self.__ITITER

    @property
    def ISLDIA(self):
        return self.__ISLDIA

    @property
    def CONVCR(self):
        return self.__CONVCR

    @property
    def ITMAX(self):
        return self.__ITMAX

    @property
    def ILUMP(self):
        return self.__ILUMP

    @property
    def NCPROJ(self):
        return self.__NCPROJ

    @property
    def NCINST(self):
        return self.__NCINST

    @property
    def NCSOUR(self):
        return self.__NCSOUR

    @property
    def NCHIST(self):
        return self.__NCHIST

    @property
    def NCREF(self):
        return self.__NCREF

    @property
    def NCCOM(self):
        return self.__NCCOM

    @property
    def NCHOST(self):
        return self.__NCHOST

    @property
    def NCCONV(self):
        return self.__NCCONV

    @property
    def NCCONT(self):
        return self.__NCCONT

    @property
    def NCDATE(self):
        return self.__NCDATE

    @property
    def FortranNamelists(self):
        return self.__FortranNamelists

    @RUNDES.setter
    def RUNDES(self, RUNDES):
        if RUNDES is None:
            self.__RUNDES = self.AdcircMesh.description.strip('\n')
        else:
            assert isinstance(RUNDES, str)
            self.__RUNDES = RUNDES.strip('\n')

    @RUNID.setter
    def RUNID(self, RUNID):
        if RUNID is None:
            today = datetime.now()
            self.__RUNID = '{}'.format(today.strftime('%Y/%m/%d'))
        else:
            assert isinstance(RUNID, str)
            self.__RUNID = RUNID.strip('\n')

    @NFOVER.setter
    def NFOVER(self, NFOVER):
        if NFOVER in [0, 1]:
            self.__NFOVER = NFOVER
        else:
            raise TypeError('NFOVER must be 0 or 1')

    @WarnElev.setter
    def WarnElev(self, WarnElev):
        if WarnElev is not None:
            self.__WarnElev = float(WarnElev)
        else:
            self.__WarnElev = None

    @iWarnElevDump.setter
    def iWarnElevDump(self, iWarnElevDump):
        if iWarnElevDump is not None:
            iWarnElevDump = int(iWarnElevDump)
            if iWarnElevDump not in [0, 1]:
                raise TypeError('iWarnElevDump must be 0 or 1')
            self.__iWarnElevDump = int(iWarnElevDump)
        else:
            if self.WarnElev is not None:
                raise RuntimeError('Must set iWarnElevDump if WarnElev is not '
                                   + 'None')

    @WarnElevDumpLimit.setter
    def WarnElevDumpLimit(self, WarnElevDumpLimit):
        if WarnElevDumpLimit is not None:
            assert isinstance(WarnElevDumpLimit, int)
            assert WarnElevDumpLimit > 0
            self.__WarnElevDumpLimit = WarnElevDumpLimit
        else:
            if self.WarnElev is not None:
                raise RuntimeError('Must set WarnElevDumpLimit if WarnElev is '
                                   + 'not None')

    @ErrorElev.setter
    def ErrorElev(self, ErrorElev):
        if ErrorElev is not None:
            self.__ErrorElev = float(ErrorElev)
        else:
            if self.WarnElev is not None:
                raise RuntimeError('Must set iWarnElevDump if WarnElev is not '
                                   + 'None')

    @NABOUT.setter
    def NABOUT(self, NABOUT):
        assert isinstance(NABOUT, int)
        assert NABOUT in [-1, 0, 1, 2, 3]
        self.__NABOUT = NABOUT

    @NSCREEN.setter
    def NSCREEN(self, NSCREEN):
        assert isinstance(NSCREEN, int)
        self.__NSCREEN = NSCREEN

    @ICS.setter
    def ICS(self, ICS):
        if ICS is None:
            if self.AdcircMesh.SpatialReference.IsProjected():
                self.__ICS = 1
            else:
                self.__ICS = 2
        else:
            ICS = int(ICS)
            assert ICS in [1, 2]
            self.__ICS = ICS

    @IM.setter
    def IM(self, IM):
        assert isinstance(IM, int)
        if IM in [0, 111112]:
            self.__IM = IM
        elif IM in [1, 2, 611112]:
            raise NotImplementedError('3D runs not yet supported.')
        else:
            raise TypeError('IM must be 0, 1, 21, 111112 or 611112.')

    @IDEN.setter
    def IDEN(self, IDEN):
        if IDEN is not None:
            raise NotImplementedError('3D runs not yet supported.')

    @NOLIBF.setter
    def NOLIBF(self, NOLIBF):
        if NOLIBF is None:
            if 'primitive_weighting_in_continuity_equation' \
              in self.AdcircMesh.NodalAttributes:
                NOLIBF = 1
            else:
                NOLIBF = 0
        else:
            NOLIBF = int(NOLIBF)
        assert NOLIBF in [0, 1, 2]
        self.__NOLIBF = NOLIBF

    @NOLIFA.setter
    def NOLIFA(self, NOLIFA):
        NOLIFA = int(NOLIFA)
        assert NOLIFA in [0, 1, 2]
        self.__NOLIFA = NOLIFA

    @NOLICA.setter
    def NOLICA(self, NOLICA):
        NOLICA = int(NOLICA)
        assert NOLICA in [0, 1]
        self.__NOLICA = NOLICA

    @NOLICAT.setter
    def NOLICAT(self, NOLICAT):
        NOLICAT = int(NOLICAT)
        assert NOLICAT in [0, 1]
        self.__NOLICAT = NOLICAT

    @NCOR.setter
    def NCOR(self, NCOR):
        if NCOR is None:
            if self.AdcircMesh.SpatialReference.IsGeographic():
                self.__NCOR = 1
            else:
                self.__NCOR = 0
        else:
            assert NCOR in [0, 1]
            self.__NCOR = NCOR

    @NTIP.setter
    def NTIP(self, NTIP):
        NTIP = int(NTIP)
        assert NTIP in [0, 1, 2]
        self.__NTIP = NTIP

    @G.setter
    def G(self, G):
        self.__G = float(G)

    @DTDP.setter
    def DTDP(self, DTDP):
        DTDP = float(DTDP)
        DTDP = np.abs(DTDP)
        assert DTDP != 0.
        self.__DTDP = DTDP

    @STATIM.setter
    def STATIM(self, STATIM):
        if STATIM is None:
            self.__STATIM = 0.
        else:
            self.__STATIM = float(STATIM)

    @REFTIM.setter
    def REFTIM(self, REFTIM):
        if REFTIM is None:
            self.__REFTIM = 0.
        else:
            self.__REFTIM = float(REFTIM)

    @DRAMP.setter
    def DRAMP(self, DRAMP):
        if DRAMP is None:
            if self.TidalForcing:
                dt = self.start_date - self.TidalForcing.forcing_start_date
                DRAMP = ((2/3)*dt.total_seconds())/(60*60*24)
        else:
            DRAMP = float(DRAMP)
        self.__DRAMP = DRAMP

    @DRAMPExtFlux.setter
    def DRAMPExtFlux(self, DRAMPExtFlux):
        self.__DRAMPExtFlux = float(DRAMPExtFlux)

    @FluxSettlingTime.setter
    def FluxSettlingTime(self, FluxSettlingTime):
        self.__FluxSettlingTime = float(FluxSettlingTime)

    @DRAMPIntFlux.setter
    def DRAMPIntFlux(self, DRAMPIntFlux):
        self.__DRAMPIntFlux = float(DRAMPIntFlux)

    @DRAMPElev.setter
    def DRAMPElev(self, DRAMPElev):
        self.__DRAMPElev = float(DRAMPElev)

    @DRAMPTip.setter
    def DRAMPTip(self, DRAMPTip):
        self.__DRAMPTip = float(DRAMPTip)

    @DRAMPMete.setter
    def DRAMPMete(self, DRAMPMete):
        self.__DRAMPMete = float(DRAMPMete)

    @DRAMPWRad.setter
    def DRAMPWRad(self, DRAMPWRad):
        self.__DRAMPWRad = float(DRAMPWRad)

    @DUnRampMete.setter
    def DUnRampMete(self, DUnRampMete):
        if DUnRampMete is None:
            DUnRampMete = self.DRAMP
        self.__DUnRampMete = float(DUnRampMete)

    @A00.setter
    def A00(self, A00):
        self.__A00 = float(A00)

    @B00.setter
    def B00(self, B00):
        self.__B00 = float(B00)

    @C00.setter
    def C00(self, C00):
        self.__C00 = float(C00)

    @H0.setter
    def H0(self, H0):
        self.__H0 = float(H0)

    @NODEDRYMIN.setter
    def NODEDRYMIN(self, NODEDRYMIN):
        self.__NODEDRYMIN = int(NODEDRYMIN)

    @NODEWETRMP.setter
    def NODEWETRMP(self, NODEWETRMP):
        self.__NODEWETRMP = int(NODEWETRMP)

    @VELMIN.setter
    def VELMIN(self, VELMIN):
        self.__VELMIN = float(VELMIN)

    @SLAM0.setter
    def SLAM0(self, SLAM0):
        if SLAM0 is None:
            SLAM0 = np.median(self.AdcircMesh.x)
        self.__SLAM0 = float(SLAM0)

    @SFEA0.setter
    def SFEA0(self, SFEA0):
        if SFEA0 is None:
            SFEA0 = np.median(self.AdcircMesh.y)
        self.__SFEA0 = float(SFEA0)

    @FFACTOR.setter
    def FFACTOR(self, FFACTOR):
        self.__FFACTOR = float(FFACTOR)

    @HBREAK.setter
    def HBREAK(self, HBREAK):
        self.__HBREAK = float(HBREAK)

    @FTHETA.setter
    def FTHETA(self, FTHETA):
        self.__FTHETA = float(FTHETA)

    @FGAMMA.setter
    def FGAMMA(self, FGAMMA):
        self.__FGAMMA = float(FGAMMA)

    @ESLM.setter
    def ESLM(self, ESLM):
        self.__ESLM = float(ESLM)

    @ESLC.setter
    def ESLC(self, ESLC):
        self.__ESLC = ESLC

    @CORI.setter
    def CORI(self, CORI):
        if CORI is None:
            if self.NCOR == 0:
                raise Exception('Must pass CORI when NCOR=0')
            else:
                CORI = 0.
        else:
            CORI = float(CORI)
        self.__CORI = CORI

    @NTIF.setter
    def NTIF(self, NTIF):
        if self.TidalForcing is None:
            NTIF = []
        else:
            if NTIF == 'all':
                NTIF = []
                for constituent in self.TidalForcing.constituents:
                    if constituent in self.TidalForcing.major8:
                        NTIF.append(constituent)
            else:
                NTIF = list(NTIF)
        self.__NTIF = NTIF

    @ANGINN.setter
    def ANGINN(self, ANGINN):
        self.__ANGINN = float(ANGINN)

    @THAS.setter
    def THAS(self, THAS):
        if THAS is not None:
            THAS = float(THAS)
            assert THAS >= 0.
        self.__THAS = THAS

    @THAF.setter
    def THAF(self, THAF):
        if THAF is not None:
            THAF = float(THAF)
            assert THAF >= 0.
        self.__THAF = THAF

    @NHAINC.setter
    def NHAINC(self, NHAINC):
        if NHAINC is not None:
            NHAINC = int(NHAINC)
            assert NHAINC >= 0
        self.__NHAINC = NHAINC

    @FMV.setter
    def FMV(self, FMV):
        FMV = float(FMV)
        assert FMV in [0., 0.1, 1.]
        self.__FMV = FMV

    @NHSTAR.setter
    def NHSTAR(self, NHSTAR):
        if NHSTAR is None:
            if self.netcdf:
                NHSTAR = 5
            else:
                NHSTAR = 1
        else:
            NHSTAR = int(NHSTAR)
            assert NHSTAR in [0, 1, 2, 3, 5]
        self.__NHSINC = NHSTAR

    @NHSINC.setter
    def NHSINC(self, NHSINC):
        dt = self.start_date - self.forcing_start_date
        if NHSINC is None:
            NHSINC = int(dt.total_seconds()/self.DTDP)
        else:
            NHSINC = int(NHSINC)
            assert NHSINC <= int((dt.total_seconds()/self.DTDP))
        self.__NHSINC = NHSINC

    @ITITER.setter
    def ITITER(self, ITITER):
        ITITER = int(ITITER)
        assert ITITER in [1, -1]
        self.__ITITER = ITITER

    @ISLDIA.setter
    def ISLDIA(self, ISLDIA):
        ISLDIA = int(ISLDIA)
        assert ISLDIA in [0, 1, 2, 3, 4, 5]
        self.__ISLDIA = ISLDIA

    @CONVCR.setter
    def CONVCR(self, CONVCR):
        self.__CONVCR = float(CONVCR)

    @ITMAX.setter
    def ITMAX(self, ITMAX):
        self.__ITMAX = int(ITMAX)

    @ILUMP.setter
    def ILUMP(self, ILUMP):
        self.__ILUMP = int(ILUMP)

    @NCPROJ.setter
    def NCPROJ(self, NCPROJ):
        self.__NCPROJ = str(NCPROJ)

    @NCINST.setter
    def NCINST(self, NCINST):
        self.__NCINST = str(NCINST)

    @NCSOUR.setter
    def NCSOUR(self, NCSOUR):
        self.__NCSOUR = str(NCSOUR)

    @NCHIST.setter
    def NCHIST(self, NCHIST):
        self.__NCHIST = str(NCHIST)

    @NCREF.setter
    def NCREF(self, NCREF):
        self.__NCREF = str(NCREF)

    @NCCOM.setter
    def NCCOM(self, NCCOM):
        self.__NCCOM = str(NCCOM)

    @NCHOST.setter
    def NCHOST(self, NCHOST):
        self.__NCHOST = str(NCHOST)

    @NCCONV.setter
    def NCCONV(self, NCCONV):
        self.__NCCONV = str(NCCONV)

    @NCCONT.setter
    def NCCONT(self, NCCONT):
        self.__NCCONT = str(NCCONT)

    @NCDATE.setter
    def NCDATE(self, NCDATE):
        if NCDATE is None:
            if self.TidalForcing:
                NCDATE = self.TidalForcing.forcing_start_date
            else:
                NCDATE = self.start_date
        else:
            NCDATE = datetime.strptime(NCDATE, '%Y-%m-%d %H:%M:%s UTC')
        self.__NCDATE = NCDATE


class Fort15TestCase(unittest.TestCase):

    def test_Fort15(self):
        AdcircMesh = Model.AdcircMesh(os.getenv('ShinnecockInlet'), 4326,
                                      'LMSL')
        now = datetime.now()
        forecast = now + timedelta(days=1)
        TidalForcing = Model._TidalForcing(now, forecast)
        fort15 = _Fort15(AdcircMesh, TidalForcing=TidalForcing)
        fort15.get_fort15('coldstart')
