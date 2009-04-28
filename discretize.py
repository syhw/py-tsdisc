#!/opt/local/bin/python2.5

import sys, csv, getopt
import pylab 
import sg_filter
import arrayfns
import math
from fisher import FisherExactTest  
fish = FisherExactTest()
# pvalue, enrichment, evaluate, print_report (k, n, C, G)
# *  k number of objects in the selection having a given property
# * n size of the selection
# * C number of objects in the population having this property
# * G size of the population

indice = 0
number_cols = 10
plot = 0
smoothing = 0
interp = 0

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

def deriv(data, time):
    """ 
    Derivate data values in regard to time with 
    d(data[i][j])/dt = (data[i+1][j]-data[i][j]) / (time[i+1]-time[i])
    we lost the last value in the process (start is more important than end)
    """
    t = [] # t[time]
    print '>>> Derivation:',
    for i in range(len(time)-1):
        if i % 100 == 0:
            print '.',
        t.append( (float(data[i+1])-float(data[i])) \
                / (time[i+1]-time[i]) )
    t.append(t[len(t)-1])
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

def interpol(data, time):
    """ 
    Normalize an array over time with spline interpolation,
    [10,20,40] taken at times [0,1,3] 
    will become [10,20,30,40] at times [0,1,2,3]
    """
    return True


def find_inflex_points(list):
    """ 
    Find relevant inflexion points in a numerical array representing
    values (evaluations) of a function
    """
    return True


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

if plot:
    pylab.figure(figsize=(8,14))
    pylab.subplot(311)
    for i in range(len(compounds)): 
        pylab.plot(time, compounds[i])
        pylab.title(first_line[i+1])

if smoothing:
    ############ smoothing (Savitzky-Golay) ############
    if plot:
        pylab.subplot(312)
    coeff = sg_filter.calc_coeff(6,3)
    smoothed_compounds = []
    for i in range(len(compounds)): 
        smoothed_compounds.append(sg_filter.smooth(compounds[i], coeff))
        if plot:
            pylab.plot(time, smoothed_compounds[i])
            pylab.title("Smoothed " + first_line[i+1])

    ############ derivative ############
    if plot:
        pylab.subplot(313)
    coeff = sg_filter.calc_coeff(6,3,1)
    smoothed_dcompounds = []
    for i in range(len(compounds)): 
        smoothed_dcompounds.append(sg_filter.smooth(compounds[i], coeff))
        #smoothed_dcompounds.append(deriv(compounds[i], time))
        if plot:
            pylab.plot(time, smoothed_dcompounds[i])
            pylab.title("Smoothed derivative of " + first_line[i+1])

    if plot:
        pylab.savefig("smoothing.png")

if interp:
    itime = [time[0]] # always an important value
    (step, t) = find_min_delta(time)
    end_time = time[len(time)-1]
    while t <= end_time:
        itime.append(t)
        t = t + step
    icompounds = []
    for c in compounds:
        icompounds.append(arrayfns.interp(c, time, itime))
    compounds = icompounds
    time = itime

    if plot:
        pylab.figure(figsize=(8,14))
        pylab.subplot(311)
        for i in range(len(compounds)): 
            pylab.plot(time, compounds[i])
            pylab.title(first_line[i+1])

    if smoothing:
        ############ smoothing (Savitzky-Golay) ############
        if plot:
            pylab.subplot(312)
        coeff = sg_filter.calc_coeff(6,3)
        smoothed_compounds = []
        for i in range(len(compounds)): 
            smoothed_compounds.append(sg_filter.smooth(compounds[i], coeff))
            if plot:
                pylab.plot(time, smoothed_compounds[i])
                pylab.title("Smoothed " + first_line[i+1])

        ############ derivative ############
        if plot:
            pylab.subplot(313)
        coeff = sg_filter.calc_coeff(6,3,1)
        smoothed_dcompounds = []
        for i in range(len(compounds)): 
            smoothed_dcompounds.append(sg_filter.smooth(compounds[i], coeff))
            #smoothed_dcompounds.append(deriv(compounds[i], time))
            if plot:
                pylab.plot(time, smoothed_dcompounds[i])
                pylab.title("Smoothed derivative of " + first_line[i+1])

        if plot:
            pylab.savefig("after_interp.png")


#fish.pvalue()

