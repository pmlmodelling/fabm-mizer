import os
import netCDF4
import numpy

min_abundance = 0.1
nspinup = 4
animate = False

nc = netCDF4.Dataset('output.nc')
instance = 'fish'
resource_instance = 'resource'

class SizeStructure:
   def __init__(self,nc,prefix):
      self.nc = nc
      self.prefix = prefix

      # Count number of bins
      self.ncohort = 0
      while '%sNw%i' % (prefix,self.ncohort+1) in nc.variables: self.ncohort += 1

      # Get mass at bin centres
      self.mass = numpy.empty((self.ncohort,))
      for icohort in range(self.ncohort):
         ncbin = self.nc.variables['%sNw%i' % (prefix,icohort+1)]
         self.mass[icohort] = ncbin.particle_mass

      # Derived mass metrics
      self.log10mass = numpy.log10(self.mass)
      self.delta_log10mass = self.log10mass[1]-self.log10mass[0]
      self.log10mass_bounds = numpy.empty((self.log10mass.size+1,),dtype=float)
      self.log10mass_bounds[:-1] = self.log10mass-0.5*self.delta_log10mass
      self.log10mass_bounds[-1] = self.log10mass[-1]+0.5*self.delta_log10mass
      self.mass_bounds = 10.**self.log10mass_bounds
      self.delta_mass = numpy.diff(self.mass_bounds)

   def getData(self,name):
      data = numpy.empty((self.ncohort,ntime),dtype=float)
      for icohort in range(self.ncohort):
         ncvar = self.nc.variables['%s%s%i' % (self.prefix,name,icohort+1)]
         data[icohort,:] = ncvar[:,0,0]
      return data

   def pcolor(self,data,long_name,units=None,**kwargs):
      if isinstance(data,basestring):
         if units is None: units = self.nc.variables['%s%s1' % (self.prefix,data)].units
         data = self.getData(data)
      pc = pylab.pcolormesh(numdts_bounds,self.mass_bounds,data,**kwargs)
      cb = pylab.colorbar(pc)
      cb.set_label('%s (%s)' % (long_name,units))
      pylab.ylabel('mass (g)')
      pylab.grid(True)
      pylab.gca().xaxis_date()
      pylab.axis('tight')
      pylab.yscale('log')

pop = SizeStructure(nc,'%s_' % instance)
resource = SizeStructure(nc,'%s_' % resource_instance)

ntime = nc.variables['time'].size
spectrum = pop.getData('Nw')/pop.delta_mass[:,numpy.newaxis]
feeding_level = pop.getData('f')
resource_spectrum = resource.getData('Nw')/resource.delta_mass[:,numpy.newaxis]
repr = pop.getData('reproduction')
nctotal_repr = nc.variables['%s_total_reproduction' % instance]
total_repr = nctotal_repr[:,0,0]
nctime = nc.variables['time']
timeunits = nctime.units
time = nctime[:]

dts = netCDF4.num2date(time[:],timeunits)

import pylab
numdts = pylab.date2num(dts)
numdts_bounds = numpy.empty((numdts.size+1,))
numdts_bounds[:-1] = numdts-(numdts[1]-numdts[0])/2
numdts_bounds[-1] = numdts[-1]+(numdts[1]-numdts[0])/2

fig = pylab.figure(figsize=(8,5))
pylab.plot_date(numdts[nspinup*365:],total_repr[nspinup*365:],'-')
pylab.ylabel('total reproduction (%s)' % nctotal_repr.units)
pylab.grid(True)
pylab.savefig('total_reproduction.png',dpi=150)
pylab.close(fig)

import matplotlib.colors

fig = pylab.figure(figsize=(8,8))
pylab.subplot(2,1,1)
pop.pcolor(repr,'reproduction',nctotal_repr.units)
pylab.subplot(2,1,2)
spectrum_masked = numpy.ma.array(spectrum,mask=spectrum<=0.)
pop.pcolor(repr/spectrum_masked,'specific reproduction','g d-1 #-1',vmin=0)
pylab.savefig('reproduction.png',dpi=150)
pylab.close(fig)

fig = pylab.figure(figsize=(8,8))
pylab.subplot(2,1,1)
pylab.loglog(resource.mass,resource_spectrum[:,-1],'-',label='resource')
pylab.loglog(pop.mass,spectrum[:,-1],'-',label='population')
pylab.xlabel('mass (g)')
pylab.ylabel('biomass density (-)')
pylab.grid(True)
pylab.legend()
pylab.ylim(1e-9,1e9)
pylab.xlim(1e-5,1e7)
pylab.subplot(2,1,2)
pylab.semilogx(pop.mass,feeding_level[:,-1],'-')
pylab.xlabel('mass (g)')
pylab.ylabel('feeding level (-)')
pylab.grid(True)
pylab.ylim(0.,1.)
pylab.xlim(1e-3,1e7)
pylab.yticks(numpy.linspace(0.,1.,9))
pylab.savefig('final_spectrum.png',dpi=150)
pylab.close(fig)

fig = pylab.figure(figsize=(8,5))
pop.pcolor(feeding_level,'feeding level','-')
pylab.savefig('f.png',dpi=150)
pylab.close(fig)

fig = pylab.figure(figsize=(8,5))
pop.pcolor('g','individual growth')
pylab.savefig('g.png',dpi=150)
pylab.close(fig)

fig = pylab.figure(figsize=(8,8))
pop.pcolor(spectrum,'biomass density','-',norm=matplotlib.colors.LogNorm(vmin=min_abundance))
pylab.savefig('spectrum.png',dpi=150)
pylab.close(fig)

fig = pylab.figure(figsize=(8,5))
for itime in range(nspinup*365+30,ntime,60):
   pylab.loglog(pop.mass,spectrum[:,itime],'-',label='day %i' % (itime-nspinup*365))
pylab.legend()
pylab.grid(True)
pylab.ylim(1e-2,spectrum.max())
#pylab.xlim(None,1e-1)
pylab.xlabel('mass (g)')
pylab.ylabel('abundance (#)')
pylab.savefig('timeseries.png',dpi=150)
pylab.close(fig)

bar_margin = 0
if animate:
   print 'Creating animation...'
   if not os.path.isdir('ani'): os.mkdir('ani')
   fig = pylab.figure()
   pylab.grid(True)
   pylab.xlabel('mass (g)')
   pylab.ylabel('abundance (#)')
   pylab.xscale('log')
   pylab.yscale('log')
   pylab.ylim(min_abundance,spectrum.max())
   #pylab.xlim(None,1e-1)
   #line, = pylab.loglog(mass,spectrum[:,0],'-')
   left = 10.**(pop.log10mass_bounds[:-1]+pop.delta_log10mass*bar_margin/2)
   right = 10.**(pop.log10mass_bounds[1:]-pop.delta_log10mass*bar_margin/2)
   rects = pylab.bar(left, spectrum[:,0], align='edge', width=right-left, color='grey')
   for itime in range(0,ntime,1):
      dt = dts[itime].strftime('%d %B')
      if itime%10==0: print '   still %i (%s)...' % (itime,dt)
      #line.set_ydata(spectrum[:,itime])
      for rect, h in zip(rects, spectrum[:,itime]): rect.set_height(h)
      pylab.title(dt)
      pylab.savefig('ani/still%04i.png' % itime,dpi=150)
   pylab.close(fig)

pylab.show()
