from utils import *
from scipy import signal

TH_AREA_NAME = {'LGd', 'LP', 'MGv', 'SGN', 'PoT', 'TH', 'Eth', 'POL', 'VPM',
               'LGv', 'PP', 'PIL', 'VPL', 'IGL', 'SGN', 'IntG', 'LD', 'MGm',
               'MGd', 'PF', 'RT'}

HIP_AREA_NAME = {'ProS', 'SUB', 'CA3', 'CA1', 'DG', 'HPF'}
VIS_AREA_PREFIX = 'VIS'


FILTERS = DataContainer()

b, a = signal.butter(8, 0.1)
FILTERS.flt1 = DataContainer()
FILTERS.flt1.Wn = 0.1
FILTERS.flt1.N = 8
FILTERS.flt1.b = b
FILTERS.flt1.a = a
FILTERS.flt1.btype = 'lowpass'
FILTERS.flt1.discription = 'filter for binned neural spiking counts with 0.2s bin size'

b, a = signal.butter(3, 0.05)
FILTERS.flt2 = DataContainer()
FILTERS.flt2.Wn = 0.05
FILTERS.flt2.N = 3
FILTERS.flt2.b = b
FILTERS.flt2.a = a
FILTERS.flt2.btype = 'lowpass'
FILTERS.flt2.discription = 'filter for running speed'

lfpsr = 1250
b, a = signal.butter(3, Wn =[.5/lfpsr*2, 4/lfpsr*2], btype = "bandpass")
FILTERS.flt3 = DataContainer()
FILTERS.flt3.Wn = [.5/lfpsr*2, 4/lfpsr*2]
FILTERS.flt3.N = 3
FILTERS.flt3.b = b
FILTERS.flt3.a = a
FILTERS.flt3.btype = 'bandpass'
FILTERS.flt3.discription = 'bandpass filter for generating delta power'

b, a = signal.butter(3, Wn = (.5+4)/2/np.pi/lfpsr*2, btype = "lowpass")
FILTERS.flt4 = DataContainer()
FILTERS.flt4.Wn = (.5+4)/2/np.pi/lfpsr*2
FILTERS.flt4.N = 3
FILTERS.flt4.b = b
FILTERS.flt4.a = a
FILTERS.flt4.btype = 'lowpass'
FILTERS.flt4.discription = 'lowpass filter for generating delta power'
