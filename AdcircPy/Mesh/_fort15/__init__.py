# from AdcircPy.Fort15 import _Fort15

class fort15(object):
    
    def __init__(self, **kwargs):
        self.RUNDES           = kwargs.pop("RUNDES", None)
        self.RUNID            = kwargs.pop("RUNID", None)
        self.NFOVER           = kwargs.pop("NFOVER", 0)
        self.NABOUT           = kwargs.pop("NABOUT", 1)
        self.NSCREEN          = kwargs.pop("NSCREEN", 100)
        self.IHOT             = kwargs.pop("IHOT", 0)
        self.ICS             = kwargs.pop("ICS", 0)
        self.IHOT             = kwargs.pop("IHOT", 0)
        self.ICS             = kwargs.pop("IHOT", 0)
        self.IM             = kwargs.pop("IHOT", 0)
        # self.IDEN 
        # self.NOLIBF
        # self.NOLIFA
        # self.NOLICA
        # self.NOLICAT
        self.reference_date   = kwargs.pop("reference_date", None)
        self.reference_time   = kwargs.pop("reference_time", None)
        self.RNDAY            = kwargs.pop("RNDAY", None)
        self.DRAMP            = kwargs.pop("DRAMP", None)
        self.nodal_attributes = kwargs.pop("nodal_attributes", None)
        self.constituent_list = kwargs.pop("constituent_list", None)
        self.boundary_TPXO    = kwargs.pop("boundary_TPXO", None)
        self._mode            = kwargs.pop("mode", None)
        self.eta_stations     = kwargs.pop("eta_stations", None)
        self.met_stations     = kwargs.pop("met_stations", None)
        self.vel_stations     = kwargs.pop("vel_stations", None)
        self.eta_start        = kwargs.pop("eta_start", None)
        self.eta_stop         = kwargs.pop("eta_stop", None)
        self.eta_steps        = kwargs.pop("eta_steps", None)
        self.vel_start        = kwargs.pop("vel_start", None)
        self.vel_stop         = kwargs.pop("vel_stop", None)
        self.vel_steps        = kwargs.pop("vel_steps", None)
        self.hotstart_steps   = kwargs.pop("hotstart_steps", None)
        self.output_filepath  = kwargs.pop("output_filepath", None)

    # def generate_forcing_from_TPXO(self, Mesh):
    #     return _Fort15.generate_forcing_from_TPXO(self, Mesh)

    
