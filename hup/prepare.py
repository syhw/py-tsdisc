#!/usr/bin/python
"""
To be compliant with GHUP:
    You provide one name and data is assumed to be in name+'_seq.csv'
    The names of the columns are normally the header of the csv file

__Inputs__:
The input file is supposed to be concentrations of chemical products over time,
the first column of the csv beeing time.
    -s : works with the speeds of changes of concentration
    -a : works with the accelerations of changes of concentration
    -p : plot

__Output__:
The output file is named 'name_disc.sol'
The format of this file is :
    conc(product, level, time).
    speed(product, level, time).
    acce(product, level, time).
With :
    - 'product' being taken in the first line of the csv file
    - 'level' being taken in '-k' levels (ordered, ex. 2 > 1 > 0)
    - 'time' being and integer taken in the discretized steps computed here

"""
# Copyright (C) 2009
# Author: Gabriel Synnaeve
# License: http://www.opensource.org/licenses/PythonSoftFoundation.php

import sys, os, subprocess, string, csv
try:
    import psyco
    psyco.log()
    #psyco.bind(deriv)
    psyco.profile(0.1)
except:
    print 'import psyco failed, running without'

### All the parameters that can be changed
prog = './gvhmm'
slice = '5'             # number of levels, 3:9:1 means between 3 and 9, incr 1
number_cols = '10'      # number of columns
indice = '0'            # indice of the ID (first indice)
times = '3'             # number of runs
viterbi_times = '2'     # number of viterbi walks
max_number_time_step = 10
predicate = 'conc' # conc / speed / acce // DO NOT SET
name = ""
for arg in sys.argv[1:]: # may use getopt
    if arg[0] != '-':
        name = arg
if name == "":
    sys.exit('Error: please give a filename as input')

def deriv(data, time):
    """ 
    Derivate data values in regard to time with 
    d(data[i][j])/dt = (data[i+1][j]-data[i][j]) / (time[i+1]-time[i])
    we lost the last value in the process (start is more important than end)
    """
    t = [] # t[time][compound_ind]
    print '>>> Derivation:',
    for i in range(len(time)-1):
        if i % 100 == 0:
            print '.',
        tmp = []
        for j in range(len(data[i])):
            tmp.append( str((float(data[i+1][j])-float(data[i][j])) \
                    / time[i+1]-time[i]) )
        t.append(tmp)
    t.append(t[len(time)-2])
    print 'Derivation finished'
    return t

def mean(t):
    mean = 0
    for e in t:
        mean = mean + int(e)
    mean = float(mean) / len(t)
    return int(round(mean))

def imin(t):
    min = 200000.0 # should use sys.maxint / float
    imin = 0
    for i in range(len(t)):
        if t[i] < min:
            min = t[i]
            imin = i
    return imin

### Read the CSV 
f = open(name+'_seq.csv')
# Try to determine the dialect and if it fails take [, / "]
try:
    dialect = csv.Sniffer().sniff(f.read(1024))
    f.seek(0)
    r = csv.reader(f, dialect)
except csv.Error, e:
    print 'Error while trying to determine the dialect of the CSV file'
    print "Warning : %s, take ',' as delimiter and '\"' as quotes" % (e)
    f.seek(0)
    r = csv.reader(f, delimiter=',', quotechar='"')
first_line = f.readline().rstrip('\n').split(',') 
#r.next() # assume than the first row is a comment
time = []
data = [] # /!\ data[time][compound_ind]
for row in r:
    time.append(float(row[0]))
    data.append(row[int(indice)+1:int(number_cols)])
f.close()

### Here we parse arguments and manipulate the values 
### if we want to do something like concentration --> speed --> acceleration
plot = 0
for arg in sys.argv:
    if arg == '-s':
        predicate = 'speed'
        data = deriv(data, time)
    elif arg == '-a':
        predicate = 'acce'
        data = deriv(deriv(data, time),time)
    elif arg == '-p':
        plot = 1

### Run the classification
list_args = [prog, '-f', name, '-k', slice, '-x', number_cols,\
        '-I', indice, '-n', times, '-z', viterbi_times]
print '>>> Will now run: '+string.join(list_args, ' ')
subprocess.call(list_args)

### Parse the Viterbi sequence to determine the levels
data_ind = range(len(data[0]))
print ">>> Will now parse the Viterbi's sequence"
f = open(name+'_viterbi.csv')
#ff = open(name+'_viterbi.csv')
try:
    dialect = csv.Sniffer().sniff(f.read(1024))
    f.seek(0)
    r = csv.reader(f, dialect)
except csv.Error, e:
    print 'Error while trying to determine the dialect of the CSV file'
    print "Warning : %s, take ',' as delimiter and '\"' as quotes" % (e)
    f.seek(0)
    r = csv.reader(f, delimiter=',', quotechar='"')

levels = [] # levels[compound_ind] gives the discretized levels sequence
levels_means = {} # levels_means[level] gives the mean value of this level
# Ugly hack time
for row in r:
    if '#### ID:LEVEL' == row[0]:
        r.next()
        r.next()
        cur = r.next()
        while ('#### ID:VITERBI-PROB' not in cur):
            levels_means[cur[0]] = cur[2]
            cur = r.next()
    elif '#### ID:VITERBI-SEQ-LEVEL-HORIZONTAL' == row[0]:
        r.next()
        r.next()
        for c in data_ind:
            levels.append(r.next()[3:])
            for j in range(int(viterbi_times) - 1):
                r.next()

### Determine the time sampling
"""
TODO: more fine grain
"""
#step_size = int(round(float(len(levels[0]))/12))
step_size = len(levels[0]) 
for c in data_ind:
    current = [c][0]
    step_size_c = 1
    min_step_size = int(round(float(len(levels[c]))/max_number_time_step))
    for i in range(len(levels[c])):
        if current != levels[c][i]:
            if step_size_c < step_size and step_size_c >= min_step_size:
                step_size = step_size_c
            current = levels[c][i] 
            step_size_c = 1
        else:
            step_size_c += 1
print '>>> Final time step size %d' % step_size

### Apply the time sampling to levels with MEAN aggreg.
fl = [] # fl[compound_ind][T] gives the discretized level from (V)HMM at T
for c in data_ind:
    number_steps = int(round(float(len(levels[c]))/step_size))
    tmp = []
    for i in range(number_steps):
        if step_size*(i+1) < len(levels[c]):
            tmp.append(mean(levels[c][step_size*i:step_size*(i+1)]))
        else:
            tmp.append(mean(levels[c][step_size*i:]))
    fl.append(tmp)
print '>>> Final levels:',
print fl

### Do plotting if necessary [DUMB DUMB DUMB]
if plot:
    import matplotlib.pyplot as plt
    data_values = [] # data[compound_ind][time]
    ### Convert data[i][c] to data_values[c][i] and set levels_values
    levels_values = []
    for c in data_ind:
        tmp = []
        tmp_lv = []
        for i in range(len(data)):
            tmp.append(float(data[i][c]))
            tmp_lv.append(float(levels_means[levels[c][i]]))
        data_values.append(tmp)
        levels_values.append(tmp_lv)
    for c in data_ind:
        plt.figure(c)
        plt.plot(time, data_values[c], 'r', linewidth=1.0)
        plt.plot(time, levels_values[c], 'k', linewidth=1.0)
        ax = [min(time)-0.05, max(time)+0.05, 
                min(levels_values[c])-0.5, max(levels_values[c])+0.5]
        plt.axis(ax)
    plt.show()

### Output the symbolic model
f = open(name+'_disc.sol', 'w')
try:
    for c in data_ind:
        for i in range(len(fl[c])):
            f.write(predicate + '(' + first_line[c+1].rstrip(',') \
                    + ',' + str(fl[c][i]) + ',' + str(i) + ').\n' )
    print '>>> Symbolic discretized model written in %s_disc.sol' % (name)
finally:
    f.close()



#if __name__ == "__main__":
#        main(sys.argv[1:])

