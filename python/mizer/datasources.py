from __future__ import print_function

import numpy
import netCDF4
from matplotlib import pyplot
from matplotlib.dates import date2num

class ValueProvider(object):
    def __init__(self):
        pass

    def get(self, time):
        raise NotImplementedError('ValueProvider-derived class MUST implement "get"')

    def mean(self):
        raise NotImplementedError('ValueProvider-derived class MUST implement "mean"')

class Constant(ValueProvider):
    def __init__(self, value):
        ValueProvider.__init__(self)
        self.value = value

    def get(self, time):
        return self.value

    def mean(self):
        return self.value

class TimeSeries(ValueProvider):
    def __init__(self, path, variable_name, scale_factor=1.0, **dim2index):
        ValueProvider.__init__(self)
        with netCDF4.Dataset(path) as nc:
            ncvar = nc.variables[variable_name]
            self.data = ncvar[...]*scale_factor
            self.units = ncvar.units
            self.long_name = ncvar.long_name
            for idim in range(ncvar.ndim-1, -1, -1):
                dimname = ncvar.dimensions[idim]
                if ncvar.shape[idim] == 1:
                    slc = [slice(None)]*self.data.ndim
                    slc[idim] = 0
                    self.data = self.data[slc]
                elif dimname in dim2index:
                    slc = [slice(None)]*self.data.ndim
                    slc[idim] = dim2index[dimname]
                    self.data = self.data[slc]
                elif dimname != 'time':
                    self.data = self.data.mean(axis=idim)
            nctime = nc.variables['time']
            self.times = date2num(netCDF4.num2date(nctime[:], nctime.units))
        self.plot()

    def get(self, time):
        return numpy.interp(time, self.times, self.data)

    def mean(self):
        return self.data.mean()

    def plot(self):
        fig = pyplot.figure()
        ax = fig.gca()
        ax.plot_date(self.times, self.data, '-')
        ax.grid()
        ax.set_ylabel('%s (g WM/m3)' % (self.long_name,))

def asValueProvider(value):
    if not isinstance(value, ValueProvider):
        value = Constant(value)
    return value
