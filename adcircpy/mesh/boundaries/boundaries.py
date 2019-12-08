
class Boundaries:

    def __init__(
        self,
        ocean_boundaries=None,
        land_boundaries=None,
        inner_boundaries=None,
        inflow_boundaries=None,
        outflow_boundaries=None,
        weir_boundaries=None,
        culvert_boundaries=None
    ):
        self._ocean_boundaries = ocean_boundaries
        self._land_boundaries = land_boundaries
        self._inner_boundaries = inner_boundaries
        self._inflow_boundaries = inflow_boundaries
        self._outflow_boundaries = outflow_boundaries
        self._weir_boundaries = weir_boundaries
        self._culvert_boundaries = culvert_boundaries

    @property
    def ocean_boundaries(self):
        return self.__ocean_boundaries

    @ocean_boundaries.setter
    def ocean_boundaries(self, ocean_boundaries):
   


        # self.ocean_boundaries = OceanBoundaries()
        # self.land_boundaries = LandBoundaries()
        # self.inner_boundaries = InnerBoundaries()
        # self.inflow_boundaries = InflowBoundaries()
        # self.outflow_boundaries = OutflowBoundaries()
        # self.weir_boundaries = WeirBoundaries()
        # self.culvert_boundaries = CulvertBoundaries()

    @property
    def ocean_boundaries(self):
        return self._ocean_boundaries
    
        land_boundaries=None,
        inner_boundaries=None,
        inflow_boundaries=None,
        outflow_boundaries=None,
        weir_boundaries=None,
        culvert_boundaries=None