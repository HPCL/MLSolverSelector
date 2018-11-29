

### Format of the CSV is [ matrix name, matrix nnz, solve time, extraction time, extraction_time, extraction ] 

import sys;
import csv
import matplotlib.pyplot as plt


def getData(filename):
  with open(filename) as csvfile:
    
    header = None
    data = {}
    data["nnz"] = []
    data["solve"] = []
    data["extract"] = []
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        if (header == None):
            header = row
        else :
            data["nnz"].append(int(row[1]))
            data["solve"].append(float(row[2]))
            for j in range(0,len(row[3:])):
                if len(data["extract"]) < j+1 :
                    data["extract"].append([])
                data["extract"][j].append( float(row[3+j]) )
    data["header"] = header
    return data;

COLORCOUNT=0;
COLORS=['b','g','r','k','m']

def plotExtractionTimes(data, ax) :
    COLORCOUNT = 0
    for i in range(0, len(data["extract"])):
        ax.scatter(data["nnz"], data["extract"][i] , c=COLORS[COLORCOUNT], alpha=0.5, label=data["header"][3+i])
        COLORCOUNT+=1
    ax.set_xscale('log')
    ax.set_xlabel("Number of non zeros");
    ax.set_ylabel("Extraction Time")

def plotSolveTimes(data, ax) :
    COLORCOUNT=0
    ax.scatter(data["nnz"], data["solve"] , alpha=0.5, label="Solve Time")
    ax.set_xlabel("Number of non zeros");
    ax.set_ylabel("Solve Time")
    ax.legend()

def plotExtractionRatio(data,ax):
    COLORCOUNT = 0
    for i in range(0, len(data["extract"])):
        ratio = []
        for j in range(0,len(data["extract"][i])):
            ratio.append( data["extract"][i][j] / data["solve"][j] ) 
        ax.scatter(data["nnz"], ratio , c=COLORS[COLORCOUNT], alpha=0.5, label=data["header"][3+i])
        COLORCOUNT += 1
    ax.legend()
    ax.set_xscale('log')
    ax.set_xlabel("Number of non zeros");
    ax.set_ylabel("Extraction Time/Solve Time")
    ax.set_xscale('log')
    ax.set_yscale('log')



plt.figure()
extract_ax = plt.gca()
plt.figure()
solve_ax = plt.gca()
plt.figure()
ratio_ax = plt.gca()

for ii in range(1, len(sys.argv)):
    data = getData(sys.argv[ii])
    plotExtractionTimes(data, extract_ax)
    plotSolveTimes(data,solve_ax)
    plotExtractionRatio(data,ratio_ax)

plt.show()



