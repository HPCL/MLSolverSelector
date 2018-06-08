

import numpy as np
import os
import csv
import hashlib
import matplotlib.pyplot as plt

# Best fit should really be "worst-best" fit. 
# Determine the linear line of best fit for the data. Should prob go cubic. 
def best_fit(X, Y):
    xbar = np.nansum(X)/len(X)
    ybar = np.nansum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = np.nansum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = np.nansum([xi**2 for xi in X]) - n * xbar**2
    b = numer / denum
    a = ybar - b * xbar

    return np.unique(X) , [ a + b * xi for xi in np.unique(X) ] 


### Extract column names :
def extract_column_names(filename):
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',') 
        header = next(reader)
        ret = [ i.strip() for i in header ] 
    return ret    

### Extract the solver information from a file. Basically this gets a 
## map from the matrix to the row that is labelled index. 
def extract_solve_map(filename, index):
    fmap = {}
    ind = 0
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',') 
        for n,row in enumerate(reader):
            if n == 0 : 
                for nn,i in enumerate(row): 
                    if i.strip() == index: 
                        ind = nn
            
            if n > 0 :
              matrix_name = os.path.splitext(os.path.basename(row[0]))[0]
              fmap[matrix_name] = float(row[ind])
    return fmap

def extract_features( filename ) :
    columns = extract_column_names(filename)
    f = {}
    for i in columns[1:]:
        f[i] = extract_solve_map(filename,i) 
    return f 


def calculate_error( filename, filename1, error_tol, plot=[] ) :
    exact_map = extract_features(filename)
    exact_map1 = extract_features(filename1);
    
    ret = {}

    count  = 0
    for i in exact_map.keys() : # i == feature 
        ret[i] = {}
        count = 0
        ret[i]['nnz'] = []
        ret[i]['error'] = []
        if i in exact_map1:  ### make sure feature in both maps. 
            c = 0
            for j in [value for value in exact_map[i] if value in exact_map1[i]] : ## matrix in both feature sets 
                c+=1
                e = abs( exact_map[i][j] - exact_map1[i][j] )/ max( 1e-16, exact_map[i][j] )
                ret[i]["nnz"].append( exact_map["nnz"][j] ) 
                ret[i]['error'].append( e ) 
                ret[i]['exact'] = exact_map[i][j] 
                ret[i]['calc'] = exact_map1[i][j] 
        ret[i]['matrices'] = c
    return ret

def plot_error_samples( exact, d, samples, plots, bestfit=True , save=False) :
    edges = [ samples[i] for i in range(0, len(samples),2)]
    inter = [ samples[i] for i in range(1, len(samples),2)]
    ax = {}
    for i in plots: 
        figure = plt.figure();
        ax[i] = plt.gca();
    
    n = 0
    colors = ["b","g","r","m","k", "b","g","r","m","k"]

    for i in range(0,len(edges)):        
        ret = calculate_error( exact, d + edges[i] + "_" + inter[i] , -1, plots );
        for j in plots:
            label = "S(" + edges[i] + "," + inter[i] + ")"
            ax[j].plot(ret[j]["nnz"] , ret[j]["error"], 'o', alpha=0.25, markeredgecolor='none', color=colors[n], label=label)
            if bestfit : 
                a,b = best_fit( ret[j]['nnz'], ret[j]['error']) 
                ax[j].plot( a,b, color = colors[n] ) 
            ax[j].set_yscale('log')
            ax[j].set_xscale('log')
            ax[j].set_xlabel('Number of non-zeros')
            ax[j].set_ylabel('Normalized error (calc - exact)/exact')
            ax[j].set_title(j)
            ax[j].legend()
        n += 1 
    if save:
        save_all_figures('Figures/error')
    plt.show()
    
def plot_extraction_time( solvefile, exactfile, samples , d, bestfit=True, save=False ) :
    edges = [ samples[i] for i in range(0, len(samples),2)]
    inter = [ samples[i] for i in range(1, len(samples),2)]

    ax = {}
    figure = plt.figure()
    ax["extract"] = plt.gca();
    figure = plt.figure()
    ax["ratio"] = plt.gca();
    figure = plt.figure()
    ax["solve"] = plt.gca();

    matrix_to_nnz = extract_solve_map(exactfile, "nnz")
    matrix_to_solve = extract_solve_map(solvefile, "Solve with GMRES+ILU")
    labels = []
    colors = ['b','r','g','m','k'] 
  
    for i in range(0,len(edges)):
        solve_time = []
        extract_time = []
        ratio = []
        nnz = []
        
        filename = d + edges[i] + "_" + inter[i] + "_timing"
        matrix_to_extract = extract_solve_map(filename, "Extract"); # extract column 4 from timing info( extract time ) 
        for matrix in matrix_to_extract:
            if matrix in matrix_to_nnz.keys() and matrix in matrix_to_solve.keys() :
                solve_time.append(matrix_to_solve[matrix])
                extract_time.append(matrix_to_extract[matrix])
                ratio.append( extract_time[-1]/solve_time[-1] )
                nnz.append( matrix_to_nnz[matrix] ) 
        label = "S(" + edges[i] + "," + inter[i] + ")"
        ax["extract"].plot( nnz, extract_time,  'o', alpha=0.25, markeredgecolor='none', color=colors[i], label=label )
        ax["ratio"].plot( nnz, ratio,  'o', alpha=0.25, markeredgecolor='none', color=colors[i], label=label )
        if bestfit:
            ae,be = best_fit( nnz, extract_time )
            ar,br = best_fit( nnz, ratio ) 
            ax["ratio"].plot( ar,br, color = colors[i] ) 
            ax["extract"].plot( ae, be, color = colors[i] )
    for j in ["ratio","extract"] :
        ax[j].set_yscale('log')
        ax[j].set_xscale('log')
        ax[j].set_xlabel('Number of non-zeros')
        ax[j].legend()
    ax["ratio"].set_ylabel('Extract time / Solve time' )
    ax["extract"].set_ylabel('Extract time ( micro seconds ) ' )
    
    if save:
        save_all_figures('Figures/extraction')
    plt.show()
    return ax


def compare_feature_sets( exactfile, samples, feature_sets, bestfit=True, labels = [], pri=False, save=False ):
    edges = [ samples[i] for i in range(0, len(samples),2)]
    inter = [ samples[i] for i in range(1, len(samples),2)]
     
    matrix_to_nnz = extract_solve_map(exactfile, "nnz")
    colors = ['b','r','g','m','k'] 
    if not len(labels) : labels = colors 
    ax = {}
    for i in range(0,len(edges)):
        plt.figure()
        ax[i] = plt.gca()
        ax[i].set_title('S('+ edges[i] + "," + inter[i] + ")" )
        ax[i].set_xlabel('Number of non-zeros')
        ax[i].set_ylabel('Extraction Time')
        for nnn,j in enumerate(feature_sets):
            extract_time = []
            nnz = []
            filename = j + edges[i] + "_" + inter[i] 
            matrix_to_extract = extract_solve_map(filename + "_timing", "Extract"); # extract column 4 from timing info( extract time ) 
            for matrix in matrix_to_extract.keys():
                if matrix in matrix_to_nnz.keys() :
                    extract_time.append(matrix_to_extract[matrix])
                    nnz.append( matrix_to_nnz[matrix] ) 
            ax[i].plot( nnz, extract_time,  'o', alpha=0.25, markeredgecolor='none', color=colors[nnn], label=labels[nnn])
            if bestfit:
                a,b = best_fit( nnz, extract_time) 
                ax[i].plot( a,b, color = colors[nnn] ) 
        ax[i].set_xscale('log')
        ax[i].set_yscale('log')
        ax[i].legend()
    if save:
        save_all_figures("Figures/comparrison_")
    plt.show()

def save_all_figures(d ):

    import matplotlib._pylab_helpers

    figures=[manager.canvas.figure
         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
    for i, figure in enumerate(figures):
        figure.savefig(d + '_%d.png' % i)

def show_all_figures():
    plt.show()

### RUN this thing
import sys

if len( sys.argv ) == 2 :
    #hardcoded stuff
    pass

else:

    # exact is the filename of the file holding the exact features 
    # estim is the filename of the file holding estimated errors 
    # error _tol is tolerance 
    # dirc is preix for sample such that dirc_i_j is a file holding features where i,j are 
    # given in the samples. 
    # plots are features to plot. 
    # feature_sets are files such that dirc_feature_set_i_j is a file 

    # python process_data.py Verify exact estim error_tol 
    if sys.argv[1] == "Verify":
        exact = sys.argv[2] 
        estim = sys.argv[3]
        error_tol = float(sys.argv[4])
        ret = calculate_error(exact, estim, error_tol ) 

        for i in ret.keys():
            
              print "%f %f %f  %s" % ( ret[i]['exact'] - ret[i]['calc'], ret[i]['exact'], ret[i]['calc'], i )

    #python process_data.py Error exact dirc/features_ 2 10 10 10 100 nnz AvgDiagonalDistance ...
    elif sys.argv[1] == "Error":
        exact = sys.argv[2]     
        dirc = sys.argv[3]


        num_samples = int(sys.argv[4])
        samples = [ sys.argv[i+5] for i in range(0,2*num_samples) ] 
        plots = sys.argv[ 5 + 2*num_samples : ] 
        
        plot_error_samples( exact, dirc, samples, plots , save=True) 
        
    #python process_data.py Extraction exact solver dir 2 10 10 10 100 
    elif sys.argv[1] == "Extraction":
        exact = sys.argv[2]     
        solver = sys.argv[3]
        dirc = sys.argv[4]
        num_samples = int(sys.argv[5])
        samples = [ sys.argv[i+6] for i in range(0,2*num_samples) ] 
        plot_extraction_time( solver, exact, samples, dirc , save=True) 
    #python process_data.py Reduced exact dirc samples SETS 
    elif sys.argv[1] == "Compare":
        
        #python process_data.py Compare $EXACT $WORKING_DIR 3 full RS1 RS2 2 10 10 10 100 
        exact = sys.argv[2]
        main = sys.argv[3]
        nsets = int(sys.argv[4])
        dircs = [ main + i for i in sys.argv[5:5+nsets] ]
        num_samples = int(sys.argv[5+nsets]) 
        samples = sys.argv[6+nsets:]
        compare_feature_sets( exact, samples, dircs, labels = [ "Full", "RS1", "RS2" ], save=True  ) 
    elif sys.argv[1] == "All"
        # This is a hardcoded example plotting everything
        #
        exact = "WORKING/features_all_interior"
        samples = [ 10 10 10 20 10 100 10 500 ]
        fsets = ["full", 'RS1', 'RS2' ] 
        fset_dircs = [ "WORKING/full/results/features_", "WORKING/RS1/results/features_", "WORKING/RS2/results/features_" ]
        compare_feature_sets( exact, samples, fset_dircs, labels =['Full','RS1','RS2'], save=True )



