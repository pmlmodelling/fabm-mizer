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
    def __init__(self, path, variable_name, scale_factor=1.0, weights=None, plot=False, **dim2index):
        ValueProvider.__init__(self)
        with netCDF4.Dataset(path) as nc:
            if variable_name in nc.variables:
                ncvar = nc.variables[variable_name]
                self.data = ncvar[...]*scale_factor
                self.long_name = ncvar.long_name
                self.units = ncvar.units
                dimensions = ncvar.dimensions
            else:
                namespace = {}
                dimensions = None
                for name, ncvar in nc.variables.items():
                    if name in variable_name:
                        namespace[name] = ncvar[...]
                        assert dimensions is None or dimensions == ncvar.dimensions
                        dimensions = ncvar.dimensions
                self.data = eval(variable_name, namespace)*scale_factor
                self.long_name = variable_name
                self.units = 'unknown'
            if scale_factor != 1.0:
                self.units = '%s*%s' % (scale_factor, self.units)
            if weights is not None:
                ncweights = nc.variables[weights]
                weight_values = ncweights[...]
                self.data *= weight_values
                self.units = '%s*%s' % (self.units, ncweights.units)
            for idim in range(len(dimensions)-1, -1, -1):
                dimname = dimensions[idim]
                if self.data.shape[idim] == 1:
                    slc = [slice(None)]*self.data.ndim
                    slc[idim] = 0
                    self.data = self.data[slc]
                    if weights is not None:
                        weight_values = weight_values[slc]
                elif dimname in dim2index:
                    if dim2index[dimname] == 'mean':
                        if weights is not None:
                            self.data = self.data.sum(axis=idim)/weight_values.sum(axis=idim)
                        else:
                            self.data = self.data.mean(axis=idim)
                    elif dim2index[dimname] == 'sum':
                        self.data = self.data.sum(axis=idim)
                    else:
                        slc = [slice(None)]*self.data.ndim
                        slc[idim] = dim2index[dimname]
                        self.data = self.data[slc]
                elif dimname != 'time':
                    assert False, 'No index (or "sum", "mean") provided for dimension %s' % dimname
            nctime = nc.variables['time']
            self.times = date2num(netCDF4.num2date(nctime[:], nctime.units))
        if plot:
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
        ax.set_ylabel('%s (%s)' % (self.long_name, self.units))

def asValueProvider(value):
    if not isinstance(value, ValueProvider):
        value = Constant(value)
    return value
