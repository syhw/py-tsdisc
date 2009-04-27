#!/opt/local/bin/python2.5

import sys, csv, getopt
import pylab 
import sg_filter
from fisher import FisherExactTest  
fish = FisherExactTest()
# pvalue, enrichment, evaluate, print_report (k, n, C, G)
# *  k number of objects in the selection having a given property
# * n size of the selection
# * C number of objects in the population having this property
# * G size of the population

pylab.figure(figsize=(8,14))
indice = 0
number_cols = 10

def Usage():
    print "./discretize.py csv_file [-d][-h][-g]"

inputname = sys.argv.pop(1)
try:    
    opts, args = getopt.getopt(sys.argv[1:],'dhg')
except getopt.GetoptError:
    Usage()
    sys.exit(2)
debug = 0
graph = 0
for o, a in opts:
    if o == '-h':
        Usage()
        sys.exit(0)
    elif o == '-d':
        debug = 1
    elif g == '-g':
        graph = 1
    else:
        assert False, "unhandled option"
#print opts
#print args

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

pylab.subplot(211)
for i in range(len(compounds)): 
    pylab.plot(time, compounds[i])
    pylab.title(first_line[i+1])

pylab.subplot(212)
coeff = sg_filter.calc_coeff(10,2)
smoothed_compounds = []
for i in range(len(compounds)): 
    smoothed_compounds.append(sg_filter.smooth(compounds[i], coeff))
    pylab.plot(time, smoothed_compounds[i])
    pylab.title("Smoothed " + first_line[i+1])
pylab.savefig("test.png")
#fish.pvalue()

