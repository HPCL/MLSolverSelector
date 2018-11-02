

### Format of the CSV is [ matrix name, matrix nnz, solve time, extraction time ] 

import sys;
import csv

filename = sys.argv[1];
header = None
data = {}
with open(filename) as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        if (header == None):
            header = row
        else :
            data[row[0]] = [ int(row[1]), float(row[2]), float(row[3]) ] ; 

import matplotlib.pyplot as plt
x = []
y = []
z = []
w = []
for k in data.keys() :
    x.append( data[k][0] )
    y.append( data[k][1] )
    z.append( data[k][2] )
    w.append( data[k][2]/data[k][1] )

f = plt.figure();
plt.scatter(x,y, alpha =0.5);
ax = plt.gca()
ax.set_xscale('log')

plt.xlabel("Number of NonZeros");
plt.ylabel("Solve Time (GMRES + ILU)");

plt.figure()
plt.scatter(x,z, alpha=0.5);
ax = plt.gca()
ax.set_xscale('log')
plt.xlabel("Number of NonZeros");
plt.ylabel("Extraction Time");

plt.figure()
plt.scatter(x,w, alpha=0.5);
ax = plt.gca()
ax.set_xscale('log')
plt.xlabel("Number of NonZeros");
plt.ylabel("Extraction Time / Solve Time ");

plt.show();




