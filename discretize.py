#!/opt/local/bin/python2.5
"""
Always reads from stdin, and print output on stdout

Default behavior:
    Read data from the provided CSV file, first column is time, 
    following ones are measured values.

-d: debug
-h: help
-p: plot
-s: activate smoothing
-i: activate interpolation

"python discretize.py data.csv"
"""
# Copyright (C) 2009 
# Author: Gabriel Synnaeve 
# License: http://www.opensource.org/licenses/PythonSoftFoundation.php

import sys, csv, getopt, warnings
import pylab 
import numpy
import sg_filter
import arrayfns
import math
from fisher import FisherExactTest  
fish = FisherExactTest()
# pvalue, enrichment, evaluate, print_report (k, n, C, G)
# * k number of objects in the selection having a given property
# * n size of the selection
# * C number of objects in the population having this property
# * G size of the population

indice = 0
number_cols = 10
plot = 0
smoothing = 0
interp = 0
first_line = ''

def deprecated(replacement):
    def new_func(func):
        func.__doc__ = "!!! Deprecated !!!" + func.__doc__
        warnings.warn(\
                "This function is deprecated, use "+replacement+" instead",\
                DeprecationWarning, 2)
        return func
    return new_func

def Usage():
    print "./discretize.py csv_file [-d][-h][-p][-s][-i]"

try:    
    inputname = sys.argv.pop(1)
    opts, args = getopt.getopt(sys.argv[1:],'dhpsi')
except getopt.GetoptError:
    Usage()
    sys.exit(2)
debug = 0
for o, a in opts:
    if o == '-h':
        Usage()
        sys.exit(0)
    elif o == '-d':
        debug = 1
    elif o == '-p':
        plot = 1
    elif o == '-s':
        smoothing = 1
    elif o == '-i':
        interp = 1 
    else:
        assert False, "unhandled option"
#print opts
#print args

def savitzky_golay(data, kernel = 11, order = 4):
    """
        applies a Savitzky-Golay filter
        input parameters:
        - data => data as a 1D numpy array
        - kernel => a positiv integer > 2*order giving the kernel size
        - order => order of the polynomal
        returns smoothed data as a numpy array

        invoke like:
        smoothed = savitzky_golay(<rough>, [kernel = value], [order = value]
    """
    try:
        kernel = abs(int(kernel))
        order = abs(int(order))
    except ValueError, msg:
        raise ValueError("kernel and order have to be of type int (floats will be converted).")
    if kernel % 2 != 1 or kernel < 1:
        raise TypeError("kernel size must be a positive odd number, was: %d" % kernel)
    if kernel < order + 2:
        raise TypeError("kernel is to small for the polynomals\nshould be > order + 2")

    # a second order polynomal has 3 coefficients
    order_range = range(order+1)
    half_window = (kernel -1) // 2
    b = numpy.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    # since we don't want the derivative, else choose [1] or [2], respectively
    m = numpy.linalg.pinv(b).A[0]
    window_size = len(m)
    half_window = (window_size-1) // 2

    # precompute the offset values for better performance
    offsets = range(-half_window, half_window+1)
    offset_data = zip(offsets, m)

    smooth_data = list()

    # temporary data, with padded zeros (since we want the same length after smoothing)
    firstval = data[0]
    lastval = data[len(data)-1]
    leftpad=numpy.zeros(half_window)+2*firstval
    rightpad=numpy.zeros(half_window)+2*lastval
    leftchunk=data[1:1+half_window]
    leftpad=leftpad-leftchunk[::-1]
    rightchunk=data[len(data)-half_window-1:len(data)-1]
    rightpad=rightpad-rightchunk[::-1]
    data = numpy.concatenate((leftpad, data))
    data = numpy.concatenate((data, rightpad))
    # data = numpy.concatenate((numpy.zeros(half_window)+firstval, data, numpy.zeros(half_window)+lastval))
    for i in range(half_window, len(data) - half_window):
        value = 0.0
        for offset, weight in offset_data:
            value += weight * data[i + offset]
        smooth_data.append(value)
    return numpy.array(smooth_data)

def plot_figs(name, time, **kwargs):
    """
    plot the figure containing all the listed kwargs, 
    values as data lists, keys as titles
    save under the given 'name'
    """
    pylab.figure(figsize=(8,4*len(kwargs)))
    num = 1
    for (k,l) in kwargs.iteritems():
        code = 100*len(kwargs) + 10 + num
        num = num + 1
        pylab.subplot(code)
        for i in range(len(l)): 
            if len(time) == len(l[i]):
                pylab.plot(time, l[i])
            else:
                div = len(time) / len(l[i]) + 1
                print div
                print len(time[::div])
                print len(l[i])
                pylab.plot(time[::div], l[i])
            pylab.title(k)
    pylab.savefig(name)

def deriv(data, time):
    """ 
    Derivate data values in regard to time with 
    d(data[i][j])/dt = (data[i+1][j]-data[i][j]) / (time[i+1]-time[i])
    we lost the last value in the process (start is more important than end)
    """
    t = [] # t[time]
    if debug:
        print '>>> Derivation:',
    for i in range(len(time)-1):
        if debug:
            if i % 100 == 0:
                print '.',
        t.append( (float(data[i+1])-float(data[i])) \
                / (float(time[i+1])-float(time[i])) )
    t.append(t[len(t)-1])
    if debug:
        print 'Derivation finished'
    return t

def find_min_delta(list):
    """ 
    Find the minimum non-nul difference between 2 consecutive elements 
    in a list and the first non-egal to the [0]-indiced element
    """
    previous_e = list[0]
    step = 100000
    found = False
    first = list[0]
    for e in list:
        if e == previous_e:
            pass
        else:
            if not found:
                found = True
                first = e
            if math.fabs(e - previous_e) < step:
                step = math.fabs(e - previous_e)
                previous_e = e 
    return (step, first)

def interpol(data, time): ### TODO spline interp
    # http://docs.scipy.org/doc/scipy/reference/interpolate.html
    # (http://www.tau.ac.il/~kineret/amit/scipy_tutorial/)
    """ 
    Normalize an array over time with spline interpolation,
    [0,1,3], [10,20,40] will become [0,1,2,3], [10,20,30,40] 
    but [ ], [ ] will become [ ], [ ]
    """
    return True

def _deriv5(h, hhx, hx, xh, xhh):
    return (-hhx - 8*hx + 8*xh + xhh)/(12*h)

def deriv5(h, list):
    demilen = 2
    constr = []
    d = []
    for i in range(demilen):
        constr.append(list[0])
    for e in list:
        constr.append(e)
    for i in range(demilen):
        constr.append(list[len(list)-1])
    for i in range(len(list)):
        d.append(_deriv5(h, constr[i], constr[i+1], \
                constr[i+3], constr[i+4]))
    return d

@deprecated("find_inflex_points")
def old_find_inflex_points(h, list):
    """ 
    Find relevant inflexion points in a numerical array representing
    values (evaluations) of a function through 5 points stencil
    """
    inflex = []
    deriv_l = deriv5(h, list)
    for i in xrange(len(list)):
        if deriv_l[i] == 0.0:
            inflex.append((i,list[i]))
    return inflex

def find_inflex_points(h, list):
    """ 
    Find relevant inflexion points in a numerical array representing
    values (evaluations) of a function
    """
    demilen = len(list) / 160
    lencontext = 2*demilen + 1
    inflex = []
    constr = []
    for i in range(demilen):
        constr.append(list[0])
    for e in list:
        constr.append(e)
    for i in range(demilen):
        constr.append(list[len(list)-1])
    context = constr[0:lencontext]
    print len(list)
    print range(len(list))
    for i in range(len(list)):
        pass
    return inflex

file = open(inputname)
r = csv.reader(file, delimiter=',', quotechar='"')
first_line = file.readline().rstrip('\n').split(',')
time = []
data = []
for row in r:
    time.append(float(row[0]))
    data.append(row[int(indice)+1:int(number_cols)])
file.close()
compounds = []
# transform data[t][c] in compounds[c][t]
for c in range(number_cols - 1):
    temp = []
    for t in xrange(len(data)):
        temp.append( float(data[t][c]) )
    t = pylab.array(temp)
    compounds.append(t)

if smoothing:
    ############ smoothing (Savitzky-Golay) ############
    coeff = sg_filter.calc_coeff(6,3)
    smoothed_compounds = []
    for i in range(len(compounds)): 
        #smoothed_compounds.append(sg_filter.smooth(compounds[i], coeff))
        smoothed_compounds.append(savitzky_golay(compounds[i]))

    ############ derivative ############
    if plot:
        pylab.subplot(313)
    coeff = sg_filter.calc_coeff(6,3,1)
    smoothed_dcompounds = []
    ttt = time[::47]
    for i in range(len(compounds)): 
        for j in range(len(compounds[i])):
            compounds[i][j] = math.log(5+compounds[i][j])
        aaa = savitzky_golay(compounds[i][::47])
        smoothed_dcompounds.append(deriv(aaa, ttt))

    if plot:
        plot_figs('smoothing.png', time, \
                compounds=compounds, \
                smoothed_compounds=smoothed_compounds, \
                smoothed_derivative_of_compounds=smoothed_dcompounds)

if interp:
    itime = [time[0]] # always an important value
    (step, t) = find_min_delta(time)
    end_time = time[len(time)-1]
    while t <= end_time:
        itime.append(t)
        t = t + step
    icompounds = []
    # linear piecewise interpolation
    for c in compounds:
        icompounds.append(arrayfns.interp(c, time, itime))
    compounds = icompounds
    time = itime

    for c in icompounds:
        old_find_inflex_points(step, c)
        pass

    if smoothing:
        ############ smoothing (Savitzky-Golay) ############
        coeff = sg_filter.calc_coeff(6,3)
        smoothed_compounds = []
        for i in range(len(compounds)): 
            smoothed_compounds.append(sg_filter.smooth(compounds[i], coeff))
        ############ derivative ############
        coeff = sg_filter.calc_coeff(6,3,1)
        smoothed_dcompounds = []
        for i in range(len(compounds)): 
            smoothed_dcompounds.append(sg_filter.smooth(compounds[i], coeff))
            #smoothed_dcompounds.append(deriv(compounds[i], time))

    if plot:
        plot_figs('after_interp.png', time, \
                compounds=compounds, \
                smoothed_compounds=smoothed_compounds, \
                smoothed_derivative_of_compounds=smoothed_dcompounds\
                )

if __name__ == '__main__':
    pass

#fish.pvalue()
