from __future__ import print_function

import os.path
import datetime
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
        return numpy.broadcast_to(self.value, numpy.shape(time))

    def mean(self):
        return self.value

class TimeSeries(ValueProvider):
    def __init__(self, path, variable_name, scale_factor=1.0, time_name='time', stop=None, minimum=None, maximum=None, allow_mask=False, expressions={}, **dim2index):
        ValueProvider.__init__(self)

        self.times = None
        print('Reading %s from %s' % (variable_name, os.path.relpath(path)))
        def getData(ncvar):
            final_dims, slc, location = [], [], []
            for dim in ncvar.dimensions:
                if dim in dim2index and isinstance(dim2index[dim], (int, slice)):
                    slc.append(dim2index[dim])
                    location.append('%s=%s' % (dim, dim2index[dim]))
                else:
                    slc.append(slice(None))
                    final_dims.append(dim)
                    location.append('%s=:' % (dim,))
            print('  %s [%s]' % (ncvar.name, ', '.join(location)))
            vardata = ncvar[slc]
            if not allow_mask:
                mask = numpy.ma.getmask(vardata)
                if mask.any():
                    first_bad_time = 'unknown' if self.times is None else num2date(self.times[mask.nonzero()[0][0]])
                    raise Exception('Variable %s in %s contains %i masked values (out of %i). Data: %s. First problem time: %s' % (path, ncvar.name, mask.sum(), mask.size, vardata, first_bad_time))
            valid = numpy.isfinite(numpy.ma.filled(vardata, 0.))
            if not valid.all():
                first_bad_time = 'unknown' if self.times is None else num2date(self.times[numpy.logical_not(valid)][0])
                raise Exception('Variable %s in %s contains %i non-finite values (NaN?). Data: %s. First problem time: %s' % (path, ncvar.name, valid.size - valid.sum(), vardata, first_bad_time))
            return final_dims, vardata

        with netCDF4.Dataset(path) as nc:
            nctime = nc.variables[time_name]
            timedim, timedata = getData(nctime)
            if stop is not None:
                istop = timedata.searchsorted(netCDF4.date2num(stop, nctime.units)) + 1
                dim2index[time_name] = slice(0, istop)
                timedata = timedata[:istop]
            self.times = date2num(netCDF4.num2date(timedata, nctime.units, only_use_cftime_datetimes=False))

            if variable_name in nc.variables:
                ncvar = nc.variables[variable_name]
                dimensions, self.data = getData(ncvar)
                self.data *= scale_factor
                self.long_name = getattr(ncvar, 'long_name', variable_name)
                self.units = getattr(ncvar, 'ncvar.units', '')
            else:
                class NcDict(object):
                    def __init__(self, nc):
                        self.nc = nc
                        self.cache = {}
                    def __getitem__(self, key):
                        if key in expressions:
                            return eval(expressions[key], {}, self)
                        if key not in self.cache:
                            ncvar = self.nc.variables[key]
                            final_dims, self.cache[key] = getData(ncvar)
                        return self.cache[key]
                    def __contains__(self, key):
                        return key in self.nc.variables

                self.data = eval(variable_name, {}, NcDict(nc))*scale_factor
                self.long_name = variable_name
                self.units = 'unknown'
                dimensions = (time_name,)
            assert self.data.shape == timedata.shape, 'Unexpected shape for %s: got %s, expected %s' % (variable_name, self.data.shape, timedata.shape)
            if scale_factor != 1.0:
                self.units = '%s*%s' % (scale_factor, self.units)
            if allow_mask:
                mask = numpy.ma.getmaskarray(self.data)
                timedata = timedata[numpy.logical_not(mask)]
                self.data = self.data.compressed()
                assert timedata.shape == self.data.shape

            self.times = date2num(netCDF4.num2date(timedata, nctime.units, only_use_cftime_datetimes=False))

        valid = numpy.isfinite(self.data)
        assert valid.all(), 'Variable %s in %s contains non-finite values (NaN?). Data: %s. First problem time: %s' % (path, variable_name, self.data, num2date(self.times[numpy.logical_not(valid)][0]))
        minval, maxval = self.data.min(), self.data.max()
        print('  Time range: %s - %s' % (num2date(self.times[0]), num2date(self.times[-1])))
        print('  Value range: %.3g - %.3g' % (minval, maxval))
        assert minimum is None or minval >= minimum, 'Minimum value %s lies below prescribed minimum of %s' % (minval, minimum)
        assert maximum is None or maxval <= maximum, 'Maximum value %s lies above prescribed maximum of %s' % (maxval, maximum)

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

class Climatology(ValueProvider):
    def __init__(self, timeseries):
        assert isinstance(timeseries, TimeSeries)
        self.times = numpy.arange(366, dtype=float)
        self.data = numpy.zeros(self.times.shape, dtype=float)
        count = numpy.zeros(self.data.shape, dtype=int)
        missing_value = timeseries.data.min() - 1.
        for iyear in range(num2date(timeseries.times[0]).year, num2date(timeseries.times[-1]).year+1):
            start = date2num(datetime.datetime(iyear, 1, 1))
            stop = date2num(datetime.datetime(iyear+1, 1, 1))
            curtime = numpy.arange(start, start + self.times[-1] + 0.1, 1.)
            istart = max(0, timeseries.times.searchsorted(start)-1)
            istop = timeseries.times.searchsorted(stop)+1
            curdata = numpy.interp(curtime, timeseries.times[istart:istop], timeseries.data[istart:istop], left=missing_value, right=missing_value)
            curcount = curdata > missing_value
            curdata[curcount == 0] = 0
            self.data += curdata
            count += curcount
        assert (count > 0).all()
        self.data /= count
        self.long_name = '%s climatology' % timeseries.long_name
        self.units = timeseries.units

    def get(self, time):
        oldshape = numpy.shape(time)
        time = numpy.reshape(time, (-1,))
        offsets = date2num([datetime.datetime(dt.year, 1, 1) for dt in num2date(time)])
        iday = time - offsets
        values = numpy.interp(iday, self.times, self.data)
        return numpy.reshape(values, oldshape)

    def mean(self):
        return self.data.mean()

    def plot(self):
        fig = pyplot.figure()
        ax = fig.gca()
        ax.plot(self.times, self.data, '-')
        ax.grid()
        ax.set_ylabel('%s (%s)' % (self.long_name, self.units))

def asValueProvider(value):
    if not isinstance(value, ValueProvider):
        value = Constant(value)
    return value
