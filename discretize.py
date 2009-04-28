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

import sys, csv, getopt
import pylab 
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

def Usage():
    print "./discretize.py csv_file [-d][-h][-p][-s][-i]"

inputname = sys.argv.pop(1)
try:    
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
            pylab.plot(time, l[i])
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

def find_inflex_points(list):
    inflex = []
    d = deriv(list, range(len(list)))
    print d
    for i in range(len(d)):
        if not d[i]:
            inflex.append(i)
    return inflex

def find_inflex_points2(list):
    """ 
    Find relevant inflexion points in a numerical array representing
    values (evaluations) of a function
    """
    inflex = []
    wl = []
    for i in range(6):
        wl.append(list[0])
    wl.append(list[:])
    for i in range(6):
        wl.append(list[len(list)-1])
    context13 = wl[0:12]
    for i in xrange(len(list)):
        context13.append(list[i+6])
        if len(context13) > 13:
            context5.pop(0)
        #if list[i+2] > list[i+1] and list[i+2] > list[i]\
        #        and list[i+2] > list[i+3] and list[i+2] > list[i+4]:
        #    inflex.append(i+2)
        #if list[i+2] < list[i+1] and list[i+2] < list[i]\
        #        and list[i+2] < list[i+3] and list[i+2] < list[i+4]:
        #    inflex.append(i+2)
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
        smoothed_compounds.append(sg_filter.smooth(compounds[i], coeff))

    ############ derivative ############
    if plot:
        pylab.subplot(313)
    coeff = sg_filter.calc_coeff(6,3,1)
    smoothed_dcompounds = []
    for i in range(len(compounds)): 
        smoothed_dcompounds.append(sg_filter.smooth(compounds[i], coeff))
        #smoothed_dcompounds.append(deriv(compounds[i], time))

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
                smoothed_derivative_of_compounds=smoothed_dcompounds)


for c in compounds:
    print c[500]
    print c[490]
#    print find_inflex_points(c)


#fish.pvalue()

