from __future__ import print_function

import os
import datetime
import netCDF4
import numpy
import pylab

class NcResult(object):
   def __init__(self, path):
      self.nc = netCDF4.Dataset(path)
   def hasVariable(self, name):
      return name in self.nc.variables
   def getUnits(self, name):
      return self.nc.variables[name].units
   def getAttribute(self, name, prop):
      return getattr(self.nc.variables[name], prop)
   def getValues(self, name):
      return self.nc.variables[name][...]
   def getTime(self):
      nctime = self.nc.variables['time']
      timeunits = nctime.units
      return netCDF4.num2date(nctime[:], timeunits, only_use_cftime_datetimes=False)

class PythonResult(object):
   def __init__(self, dts, name2data, name2var):
      self.dts = dts
      self.name2data = name2data
      self.name2var = name2var
   def hasVariable(self, name):
      return name in self.name2data
   def getUnits(self, name):
      return self.name2var[name].units
   def getAttribute(self, name, prop):
      return self.name2var[name].getRealProperty(prop)
   def getValues(self, name):
      return self.name2data[name]
   def getTime(self):
      return self.dts

class SizeStructure:
   def __init__(self, result, prefix):
      self.result = result
      self.prefix = prefix
      self.dts = self.result.getTime()
      self.ntime = len(self.dts)

      # Count number of bins
      self.ncohort = 0
      while result.hasVariable('%sNw%i' % (prefix,self.ncohort+1)):
         self.ncohort += 1

      # Get mass at bin centres
      self.mass = numpy.empty((self.ncohort,))
      for icohort in range(self.ncohort):
         self.mass[icohort] = result.getAttribute('%sNw%i' % (prefix, icohort + 1), 'particle_mass')

      # Derived mass metrics
      self.log10mass = numpy.log10(self.mass)
      self.delta_log10mass = self.log10mass[1] - self.log10mass[0]
      self.log10mass_bounds = numpy.empty((self.log10mass.size + 1,),dtype=float)
      self.log10mass_bounds[:-1] = self.log10mass-0.5*self.delta_log10mass
      self.log10mass_bounds[-1] = self.log10mass[-1]+0.5*self.delta_log10mass
      self.mass_bounds = 10.**self.log10mass_bounds
      self.delta_mass = numpy.diff(self.mass_bounds)

   def getData(self,name):
      data = numpy.empty((self.ncohort, self.ntime), dtype=float)
      for icohort in range(self.ncohort):
         data[icohort, :] = numpy.squeeze(self.result.getValues('%s%s%i' % (self.prefix, name, icohort + 1)))
      return data

   def pcolor(self, data, long_name, units=None, **kwargs):
      numdts = pylab.date2num(self.result.getTime())
      numdts_bounds = numpy.empty((numdts.size + 1,))
      numdts_bounds[:-1] = numdts-(numdts[1] - numdts[0]) / 2
      numdts_bounds[-1] = numdts[-1] + (numdts[1] - numdts[0]) / 2

      if isinstance(data,basestring):
         if units is None:
            units = self.result.getUnits('%s%s1' % (self.prefix,data))
         data = self.getData(data)
      pc = pylab.pcolormesh(numdts_bounds,self.mass_bounds,data,**kwargs)
      cb = pylab.colorbar(pc)
      cb.set_label('%s (%s)' % (long_name,units))
      pylab.ylabel('mass (g)')
      pylab.grid(True)
      pylab.gca().xaxis_date()
      pylab.axis('tight')
      pylab.yscale('log')

def plot(result, nspinup=0, instance='fish', resource_instance='resource', animate=False, spacing=30, min_density=None):
   pop = SizeStructure(result, '%s_' % instance)
   resource = SizeStructure(result, '%s_' % resource_instance)
   dts = result.getTime()
   numdts = pylab.date2num(dts)

   start = datetime.datetime(dts[0].year + nspinup, dts[0].month, dts[0].day)
   istart = numdts.searchsorted(pylab.date2num(start))

   total_repr = result.getValues('%s_total_reproduction' % instance)
   total_repr_units = result.getUnits('%s_total_reproduction' % instance)
   total_rec = result.getValues('%s_R' % instance)
   total_rec_units = result.getUnits('%s_R' % instance)
   repr = pop.getData('reproduction')
   spectrum = pop.getData('Nw') / pop.delta_mass[:, numpy.newaxis]
   feeding_level = pop.getData('f')
   resource_spectrum = resource.getData('Nw')/resource.delta_mass[:, numpy.newaxis]
   if min_density is None:
      min_density = spectrum.min()

   fig = pylab.figure(figsize=(8, 5))
   pylab.plot_date(numdts[istart:], total_repr[istart:], '-')
   pylab.ylabel('total reproduction (%s)' % total_repr_units)
   pylab.grid(True)
   pylab.savefig('total_reproduction.png', dpi=150)
   pylab.close(fig)

   fig = pylab.figure(figsize=(8, 5))
   pylab.plot_date(numdts[istart:], total_rec[istart:], '-')
   pylab.ylabel('recruitment (%s)' % total_rec_units)
   pylab.grid(True)
   pylab.savefig('recruitment.png', dpi=150)
   pylab.close(fig)

   import matplotlib.colors

   fig = pylab.figure(figsize=(8, 8))
   pylab.subplot(2, 1, 1)
   pop.pcolor(repr,'reproduction', total_repr_units)
   pylab.subplot(2, 1, 2)
   spectrum_masked = numpy.ma.array(spectrum, mask=spectrum <= 0.)
   pop.pcolor(repr / spectrum_masked, 'specific reproduction', 'g d-1 #-1', vmin=0)
   pylab.savefig('reproduction.png', dpi=150)
   pylab.close(fig)

   fig = pylab.figure(figsize=(8, 8))
   pylab.subplot(2, 1, 1)
   pylab.loglog(resource.mass, resource_spectrum[:,-1], '-', label='resource')
   pylab.loglog(pop.mass, spectrum[:,-1], '-', label='population')
   pylab.xlabel('mass (g)')
   pylab.ylabel('biomass density (-)')
   pylab.grid(True)
   pylab.legend()
   pylab.ylim(1e-9, 1e9)
   pylab.xlim(1e-5, 1e7)
   pylab.subplot(2, 1, 2)
   pylab.semilogx(pop.mass, feeding_level[:,-1], '-')
   pylab.xlabel('mass (g)')
   pylab.ylabel('feeding level (-)')
   pylab.grid(True)
   pylab.ylim(0., 1.)
   pylab.xlim(1e-3, 1e7)
   pylab.yticks(numpy.linspace(0., 1., 9))
   pylab.savefig('final_spectrum.png', dpi=150)
   pylab.close(fig)

   fig = pylab.figure(figsize=(8, 5))
   pop.pcolor(feeding_level,'feeding level', '-')
   pylab.savefig('f.png', dpi=150)
   pylab.close(fig)

   fig = pylab.figure(figsize=(8, 5))
   pop.pcolor('g', 'individual growth')
   pylab.savefig('g.png', dpi=150)
   pylab.close(fig)

   fig = pylab.figure(figsize=(8, 8))
   pop.pcolor(spectrum, 'biomass density', '-', norm=matplotlib.colors.LogNorm(vmin=min_density))
   pylab.savefig('spectrum.png', dpi=150)
   pylab.close(fig)

   fig = pylab.figure(figsize=(8, 5))
   for itime in range(istart, spectrum.shape[1], spacing):
      pylab.loglog(pop.mass, spectrum[:, itime], '-', label=dts[itime].strftime('%Y-%m-%d'))
   pylab.legend()
   pylab.grid(True)
   #pylab.ylim(1e-2, spectrum.max())
   #pylab.xlim(None,1e-1)
   pylab.xlabel('mass (g)')
   pylab.ylabel('biomass density (-)')
   pylab.savefig('timeseries.png', dpi=150)
   pylab.close(fig)

   bar_margin = 0
   if animate:
      print('Creating animation...')
      if not os.path.isdir('ani'):
         os.mkdir('ani')
      fig = pylab.figure()
      pylab.grid(True)
      pylab.xlabel('mass (g)')
      pylab.ylabel('biomass density (-)')
      pylab.xscale('log')
      pylab.yscale('log')
      pylab.ylim(min_density, spectrum.max())
      #pylab.xlim(None,1e-1)
      #line, = pylab.loglog(mass,spectrum[:,0],'-')
      left = 10.**(pop.log10mass_bounds[:-1]+pop.delta_log10mass*bar_margin/2)
      right = 10.**(pop.log10mass_bounds[1:]-pop.delta_log10mass*bar_margin/2)
      rects = pylab.bar(left, spectrum[:,0], align='edge', width=right-left, color='grey')
      for itime in range(0, len(dts), 1):
         dt = dts[itime].strftime('%d %B')
         if itime % 10 == 0:
            print('   still %i (%s)...' % (itime, dt))
         #line.set_ydata(spectrum[:,itime])
         for rect, h in zip(rects, spectrum[:,itime]): rect.set_height(h)
         pylab.title(dt)
         pylab.savefig('ani/still%04i.png' % itime, dpi=150)
      pylab.close(fig)

if __name__ == '__main__':
   plot(NcResult('output.nc'))
