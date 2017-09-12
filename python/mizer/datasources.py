from __future__ import print_function

import numpy
import netCDF4
from matplotlib import pyplot
from matplotlib.dates import date2num, num2date

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
    def __init__(self, path, variable_name, scale_factor=1.0, weights=None, plot=False, time_name='time', stop=None, **dim2index):
        ValueProvider.__init__(self)

        def getData(ncvar):
            final_dims, slc, location = [], [], []
            for dim in ncvar.dimensions:
                if dim in dim2index and isinstance(dim2index[dim], (int, slice)):
                    slc.append(dim2index[dim])
                    location.append('%s=%s' % (dim, dim2index[dim]))
                else:
                    slc.append(slice(None))
                    final_dims.append(dim)
                    location.append('%s=all' % (dim,))
            print('Reading %s from %s at %s' % (ncvar.name, path, ', '.join(location)))
            vardata = ncvar[slc]
            valid = numpy.isfinite(vardata)
            assert valid.all(), 'Variable %s in %s contains non-finite values (NaN?). Data: %s. First problem time: %s' % (path, ncvar.name, vardata, num2date(self.times[numpy.logical_not(valid)][0]))
            return final_dims, vardata

        with netCDF4.Dataset(path) as nc:
            nctime = nc.variables[time_name]
            timedim, timedata = getData(nctime)
            if stop is not None:
                istop = timedata.searchsorted(netCDF4.date2num(stop, nctime.units)) + 1
                dim2index[time_name] = slice(0, istop)
                timedata = timedata[:istop]
            self.times = date2num(netCDF4.num2date(timedata, nctime.units))

            if variable_name in nc.variables:
                ncvar = nc.variables[variable_name]
                dimensions, self.data = getData(ncvar)
                self.data *= scale_factor
                self.long_name = ncvar.long_name
                self.units = ncvar.units
            else:
                namespace = {}
                dimensions = None
                for name, ncvar in nc.variables.items():
                    if name in variable_name:
                        vardims, namespace[name] = getData(ncvar)
                        assert dimensions is None or dimensions == vardims
                        dimensions = vardims
                self.data = eval(variable_name, namespace)*scale_factor
                self.long_name = variable_name
                self.units = 'unknown'
            if scale_factor != 1.0:
                self.units = '%s*%s' % (scale_factor, self.units)
            if weights is not None:
                ncweights = nc.variables[weights]
                weight_dims, weight_values = getData(ncweights)
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
                        assert False, 'Unknown dimension indexer %s specified for %s' % (dim2index[dimname], dimname)
                elif dimname != time_name:
                    assert False, 'No index (or "sum", "mean") provided for dimension %s' % dimname
        valid = numpy.isfinite(self.data)
        assert valid.all(), 'Variable %s in %s contains non-finite values (NaN?). Data: %s. First problem time: %s' % (path, variable_name, self.data, num2date(self.times[numpy.logical_not(valid)][0]))
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
