#!/usr/bin/env python

# Extracts initial condition for online GOTM-mizer from offline result

import argparse
import netCDF4

def convert(ncpath, yamlpath, itime, scale):
    with netCDF4.Dataset(ncpath) as nc, open(yamlpath, 'w') as fout:
        ws = nc['w'][:]
        spectrum = scale * nc['spectrum'][itime, :]
        for i, (w, value) in enumerate(zip(ws, spectrum)):
            fout.write('c%i: %s   # biomass (%s * g/m2) for bin with individual mass of %.3g g\n' % (i + 1, value, scale, w))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('result')
    parser.add_argument('outfile')
    parser.add_argument('--itime', type=int, default=0)
    parser.add_argument('--scale', type=float, default=1.)
    args = parser.parse_args()
    convert(args.result, args.outfile, args.itime, scale=args.scale)