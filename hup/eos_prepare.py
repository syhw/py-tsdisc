#!/usr/bin/python
# vim: set fileencoding=utf-8
# encoding=utf-8
"""
To be compliant with HUP:
    You provide one name and data is assumed to be in name+'_seq.csv'
    The names of the columns are normally the header of the csv file

__Inputs__:
The input file is supposed to be concentrations of chemical products over time,
the first column of the csv beeing time.
/!\ Time is considered "as is" and dealt with like if equal-width
    -s : works with the speeds of changes of concentration
    -a : works with the accelerations of changes of concentration

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
prog = './vhmmd'
slice = '3'             # number of levels, 3:9:1 means between 3 and 9, incr 1
number_cols = '10'      # number of columns
indice = '0'            # indice of the ID (first indice)
times = '3'             # number of runs
viterbi_times = '3'     # number of viterbi walks
max_number_time_steps = 10
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
                    / (float(time[i+1])-float(time[i]))) )
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
r.next() # assume than the first row is a comment
time = []
data = [] # data[time][compound_ind]
for row in r:
    time.append(row[0])
    data.append(row[int(indice)+1:int(number_cols)])
f.close()

### Here we manipulate the values if we want to do something like
### concentration --> speed --> acceleration
for arg in sys.argv:
    if arg == '-s':
        predicate = 'speed'
        data = deriv(data, time)
    elif arg == '-a':
        predicate = 'acce'
        data = deriv(deriv(data, time),time)

### Write the new CSV that we will use with HUP
eosname = 'eos_'+name
f = open(eosname+'_seq.csv', 'w')
w = csv.writer(f, delimiter=',', quotechar='"')
data_ind = range(len(data[0]))
for i in data_ind:
    for j in range(len(time)):
        w.writerow([time[j],data[j][i]])
    w.writerow(['EOS'])
f.close()
print '>>> The new CSV has been written in %s_seq.csv' % (eosname)

### Run the classification
list_args = [prog, '-f', eosname, '-k', slice, '-x', '2',\
        '-I', indice, '-n', times, '-z', viterbi_times]
list_args.append('-u')
list_args.append('1')
print '>>> Will now run: '+string.join(list_args, ' ')
subprocess.call(list_args)

### Parse the Viterbi sequence to determine the levels
print ">>> Will now parse the Viterbi's sequence"
f = open(eosname+'_viterbi.csv')
try:
    dialect = csv.Sniffer().sniff(f.read(1024))
    f.seek(0)
    r = csv.reader(f, dialect)
except csv.Error, e:
    print 'Error while trying to determine the dialect of the CSV file'
    print "Warning : %s, take ',' as delimiter and '\"' as quotes" % (e)
    r = csv.reader(f, delimiter=',', quotechar='"')
levels = [] # levels[compound_ind] gives the discretized levels sequence
# Ugly hack time
#stack_ind = data_ind[::]
enc = 0
for row in r:
    if enc:
        #if int(row[0]) == stack_ind[0] and row[1] == '0':
        #    stack_ind.pop(0)
        levels.append(row[2:])
        for j in range(int(viterbi_times) - 1):
            r.next()
        #else:
        #    print "Viterbi's sequence parsing failed"
    if '#### ID:VITERBI-SEQ-HORIZONTAL' == row[0]:
        enc = 1
        r.next()
        r.next()
        
### Determine the time sampling and apply it to levels with MEAN aggreg.
"""
TODO: more fine grain
"""
#step_size = int(round(float(len(levels[0]))/12))
step_size = len(levels[0]) 
for c in data_ind:
    current = [c][0]
    step_size_c = 1
    min_step_size = int(round(float(len(levels[c]))/max_number_time_steps))
    for i in range(len(levels[c])):
        if current != levels[c][i]:
            if step_size_c < step_size and step_size_c >= min_step_size:
                step_size = step_size_c
            current = levels[c][i] 
            step_size_c = 1
        else:
            step_size_c += 1
print '>>> Final time step size %d' % step_size
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

### Reorder with level j > level i when j > i
if 'vhmm' in prog:
    param = '_hparam'
    stringgauss = 'PARAM-GAUSS-MEAN0'
else:
    param = '_param'
    stringgauss = 'PARAM-EMIT-GAUSS-MEAN0'
f = open(eosname + param + '.csv')
r = csv.reader(f, delimiter=',', quotechar='"')
# Ugly hack time #2
levelvalue = []
for row in r:
    if stringgauss in row[0]:
        enc = 1
        r.next()
        r.next()
        currow =  r.next()
        while '###' not in currow[0]  :
            levelvalue.append(float(currow[3])) # Assume order in _hparam.csv
            currow =  r.next()
        break
f.close()
maplevels = {}
count = 0
filler = max(levelvalue) + 10 # arbitrary
for i in range(len(levelvalue)):
    im = imin(levelvalue)
    maplevels[str(im)] = count
    levelvalue[im] = filler
    count += 1
time_steps = range(len(fl[0]))
print '>>> Unmapped levels:',
print fl
for c in data_ind:
    for t in time_steps:
        fl[c][t] = maplevels[str(fl[c][t])]
print '>>> Final levels:',
print fl

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

