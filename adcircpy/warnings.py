import warnings

class ADCIRCPyWarning(UserWarning):
    '''
    base adcircpy warning class
    '''
    pass

class UnsupportedFeatureWarning(ADCIRCPyWarning):
    '''
    Feature exists in ADCIRC but is not yet supported by adcircpy
    '''
    pass

class ModelSetupWarning(ADCIRCPyWarning):
    '''
    Indicates a configuration or setup issue in adcircpy inputs
    '''
    pass

def warn_adcirc(message, category=ADCIRCPyWarning, stacklevel=2):
    '''
    Prints a standardized warning message

    Parameters
    ----------
    message: str
        The warning message to display
    category: Warning subclass, optional
        The type of warning to issue (defaults to ADCIRCPyWarning)
    stacklevel : int, optional
        The stack level for the warning location
    '''
    clean_message = (message.strip()).rstrip(".") + "."
    warnings.warn(clean_message, category, stacklevel=stacklevel)
