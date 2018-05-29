from AdcircPy.Fort15 import _Fort15

class Fort15(object):
    
    def __init__(self, **kwargs):
        _Fort15.__init__(self, **kwargs)

    def generate_forcing_from_TPXO(self, Mesh):
        return _Fort15.generate_forcing_from_TPXO(self, Mesh)

    
