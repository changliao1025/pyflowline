class pycell(object):
    __metaclass__ = ABCMeta 
    pVertex_start = None
    pVertex_end = None
    dLength=0.0

    def __init__(self, pVertex_start_in, pVertex_end_in):
        self.pVertex_start = pVertex_start_in
        self.pVertex_end = pVertex_end_in
        return

    def calculate_length(self):
        dLength =0.0

        dLength = self.pVertex_start.calculate_distance( self.pVertex_end)
        self.dLength= dLength

        return dLength