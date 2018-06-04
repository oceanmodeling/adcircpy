from AdcircPy.Validation import _HighWaterMarks

class HighWaterMarks(dict):
    
    def __init__(self, **kwargs):
        dict.__init__(self, **kwargs)

    @staticmethod
    def from_csv(path):
        return _HighWaterMarks.from_csv(path)

    @staticmethod
    def from_event_id(event_id):
        _HighWaterMarks.from_event_id(event_id)

    def remove(self, type_list, **kwargs):
        """
        type_list is a list or a string. 
        Allowed strings in type_list are: 'Excellent', 'Good', 'Fair', 'Poor'
        Strings in type_list are case insensitive.
        """
        return _HighWaterMarks.remove(self, type_list, **kwargs)

    def clip_from_shapefile(self, path, **kwargs):
        return _HighWaterMarks.clip_from_shapefile(self, path, **kwargs)

    def convert_to_meters(self):
        _HighWaterMarks.convert_to_meters(self)

    def get_coordinates(self):
        return _HighWaterMarks.get_coordinates(self)

    def get_values(self):
        return _HighWaterMarks.get_values(self)
    
    def get_xyz(self):
        return _HighWaterMarks.get_xyz(self)
    
    def get_environments(self):
        return _HighWaterMarks.get_environments(self)
        
    def get_counties(self):
        return _HighWaterMarks.get_counties(self)
        
    def still_water_only(self, **kwargs):
        return _HighWaterMarks.still_water_only(self, **kwargs)

    def export_shapefile(self, path):
        _HighWaterMarks.export_shapefile(self, path)

    def make_plot(self, **kwargs):
        return _HighWaterMarks.make_plot(self, **kwargs)