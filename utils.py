import numpy as np
import random

""" Define data structure """
class StimulusSpecificData:
    def __init__(self, name):
        self.stimName = name
        
class DataContainer:
    def __init__(self):
        pass

class Ensembles:
    def __init__(self):
        pass
    
    def esmbl(self, esmblID):
        attr = 'esmbl{:d}'.format(esmblID)
        return getattr(self, attr)
    
    def setup_esmbl(self, esmblID, obj):
        attr = 'esmbl{:d}'.format(esmblID)
        setattr(self, attr, obj)


class AttrDict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def set_random_seed():
    np.random.seed(1234)
    random.seed(1234)