import sys
import os
import datetime

import numpy
import scipy.integrate

import plot_distribution

def simulate():

    # Create model (loads fabm.yaml)
    model = pyfabm.Model()

    # Configure the environment
    # Note: the set of environmental dependencies depends on the loaded biogeochemical model.
    model.findDependency('bottom_depth').value = 1.
    model.findDependency('fish/interaction_depth').value = 1.

    # Verify the model is ready to be used
    assert model.checkReady(), 'One or more model dependencies have not been fulfilled.'

    # Time derivative
    def dy(y, t0):
        model.state[:] = y
        return model.getRates()

    def get_diagnostics(y, t0):
        model.state[:] = y
        model.getRates()
        return numpy.array([variable.value for variable in model.diagnostic_variables if variable.output])

    # Time-integrate (note: FABM's internal time unit is seconds!)
    t = numpy.arange(0., 365. * 100, 5.)
    print('Simulating...')
    y = scipy.integrate.odeint(dy, model.state, t*86400)

    # Obtain diagnostics for each output time
    print('Computing diagnostics...')
    y_diag = numpy.array([get_diagnostics(y[i, :], t[i]) for i in range(y.shape[0])])

    # Convert times into dates
    start = datetime.datetime(2000, 1, 1)
    dts = [start + datetime.timedelta(days=d) for d in t]

    # Gather outputs
    name2data, name2var = {}, {}
    i = 0
    for variable in model.diagnostic_variables:
        if variable.output:
            name2data[variable.output_name] = y_diag[:, i]
            name2var[variable.output_name] = variable
            i += 1
    for i, variable in enumerate(model.state_variables):
        name2data[variable.output_name] = y[:, i]
        name2var[variable.output_name] = variable

    print('Plotting...')
    res = plot_distribution.PythonResult(dts, name2data, name2var)
    plot_distribution.plot(res, spacing=int(365/5.*5), animate=True)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pyfabmdir', action='append')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()
    for dir in args.pyfabmdir:
        sys.path.insert(0, dir)
    if args.debug:
        print('Process id = %s. Press enter to start...' % os.getpid())
        raw_input()
    import pyfabm
    simulate()
