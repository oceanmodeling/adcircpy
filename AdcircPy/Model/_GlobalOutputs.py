from datetime import timedelta


class _GlobalOutputs(object):
    def __init__(self, sampling_frequency=timedelta(0), spinup=False,
                 harmonic_analysis=False, netcdf=True):
        self._sampling_frequency = sampling_frequency
        self._spinup = spinup
        self._harmonic_analysis = harmonic_analysis
        self._netcdf = netcdf

    def __call__(self, runtype, forcing_start_date, start_date,
                 end_date, DTDP):
        return self.__get_fort15_params(runtype, forcing_start_date,
                                        start_date, end_date, DTDP)

    def _get_NHA(self, runtype):
        assert runtype in ['coldstart', 'hotstart', 'metonly']
        if runtype == 'coldstart':
            if self.spinup and self.harmonic_analysis:
                return 1
            else:
                return 0
        elif runtype == 'hotstart':
            if self.harmonic_analysis:
                return 1
            else:
                return 0
        else:
            NotImplementedError('Can we do harmonic analysis on winds only for'
                                + 'metonly runs?')
            return 0

    def __get_fort15_params(self, runtype, forcing_start_date, start_date,
                            end_date, DTDP):
        NSPOOL = int(self.sampling_frequency.total_seconds()/DTDP)
        if runtype == 'coldstart':
            if self.spinup:
                if NSPOOL > 0.:
                    if self.netcdf:
                        NOUT = -5
                    else:
                        NOUT = -1
                    TOUTST = 0.
                    dt = start_date - forcing_start_date
                    TOUTF = dt.total_seconds() / (24.*60.*60.)
                else:
                    NOUT = 0
                    TOUTST = 0.
                    TOUTF = 0.
                    NSPOOL = 0
            else:
                NOUT = 0
                TOUTST = 0.
                TOUTF = 0.
                NSPOOL = 0

        elif runtype == 'hotstart':
            if NSPOOL > 0.:
                if self.netcdf:
                    NOUT = -5
                else:
                    NOUT = -1
            else:
                NOUT = 0
            dt = start_date - forcing_start_date
            TOUTST = dt.total_seconds() / (24.*60.*60.)
            dt = end_date - forcing_start_date
            TOUTF = dt.total_seconds() / (24.*60.*60.)
        return NOUT, TOUTST, TOUTF, NSPOOL

    @property
    def sampling_frequency(self):
        return self._sampling_frequency

    @property
    def spinup(self):
        return self._spinup

    @property
    def harmonic_analysis(self):
        return self._harmonic_analysis

    @property
    def netcdf(self):
        return self._netcdf

    @property
    def NOUT(self):
        return self.__NOUT

    @property
    def _sampling_frequency(self):
        return self.__sampling_frequency

    @property
    def _spinup(self):
        return self.__spinup

    @property
    def _harmonic_analysis(self):
        return self.__harmonic_analysis

    @property
    def _netcdf(self):
        return self.__netcdf

    @_sampling_frequency.setter
    def _sampling_frequency(self, sampling_frequency):
        assert isinstance(sampling_frequency, timedelta)
        self.__sampling_frequency = sampling_frequency

    @_spinup.setter
    def _spinup(self, spinup):
        assert isinstance(spinup, bool)
        self.__spinup = spinup

    @_harmonic_analysis.setter
    def _harmonic_analysis(self, harmonic_analysis):
        assert isinstance(harmonic_analysis, bool)
        self.__harmonic_analysis = harmonic_analysis

    @_netcdf.setter
    def _netcdf(self, netcdf):
        assert isinstance(netcdf, bool)
        self.__netcdf = netcdf
