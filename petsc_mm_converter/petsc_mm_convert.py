import PetscBinaryIO as PBIO
import scipy.io
import sys, os 
import numpy as np
import csv
from subprocess import call

def get_group_map(stats_file='ssstats.csv'):
    group_map = {}
    with open( stats_file, "rb") as csvfile :
        reader = csv.reader(csvfile, delimiter=',')
        for n,row in enumerate(reader):
            if n > 2: 
                group_map[row[1]] = row[0] 
    return group_map 

def convert_mm_binary(filename, output):
    A = scipy.io.mmread(filename)
    io = PBIO.PetscBinaryIO()
    io.writeMatSciPy(open(output,'w'),A) 

def extract_tar( filename, edir="./"):
    return call(["tar", "-xf", filename, "-C" , edir ] )

def download_matrix(filename, group, f, d ) :
    dd = "https://sparse.tamu.edu/" + f + "/" + group + "/" + filename  
    return call(["wget", dd , "-O", d + filename ])

def delete_file( filename ):
    return call(["rm",filename])
def delete_dir( dire ):
    return call(["rm", "-r", dire])

format_map = { "MM" : ".mtx" , "RB" : ".rb", "matlab" : ".mat" } 

def download_extract_convert( matrix, d, extension, f ) :
    mname = matrix + ".tar.gz"
    dname = d;
    gmap = get_group_map()
    if matrix not in gmap:
        return 1 
    group = gmap[matrix] 
    r = download_matrix(mname, group, f, dname)
    if not r:
        r = extract_tar( dname + mname, d )
    print r
    if not r:
        convert_mm_binary( dname + matrix + "/" + matrix + format_map[f] , dname + matrix + extension )  
    if not r:
        r = delete_file( dname + mname )
    if not r:
        r = delete_dir(dname + matrix )
    return r 

def download_if_not_exists( matrix, d="./MM/", extension='.pbin', f="MM" ):
    import os.path 
    if not os.path.isfile(d + matrix + extension ):
        return download_extract_convert( matrix, d, extension, f )
    else :
        print "Already exists"
        return 2
def download_based_on_names( filename, direct ):
    with open( filename, "r") as csvfile :
        reader = csv.reader(csvfile, delimiter=',')
        for n,row in enumerate(reader):
            if n>0 : 
                matrix_file = os.path.basename(row[0])
                matrix_name = os.path.splitext(matrix_file)[0]
                r = download_if_not_exists(matrix_name, d=direct)
                if r == 0 : print "Successfully downloaded " , matrix_name 
                elif r == 1 : print "download failed ", matrix_name
                else : print "Already downloaded" , matrix_name
                
download_based_on_names(sys.argv[1], sys.argv[2])



