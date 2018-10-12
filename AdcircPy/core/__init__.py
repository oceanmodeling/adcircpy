"""
Function for aliasing kwargs used in _AdcircRun.py
Reference:
https://stackoverflow.com/questions/29374425/aliasing-in-the-names-of-function-arguments-double-naming
"""
import functools
def alias(aliases):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*argv, **kwargs):
            for name, alias in aliases.items():
                if name not in kwargs and alias in kwargs:
                    kwargs[name] = kwargs[alias]
            return func(*argv, **kwargs)
        return wrapper
    return decorator