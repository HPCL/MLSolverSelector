import PetscBinaryIO as PBIO
import scipy.io
import sys, os 

filename = sys.argv[1] 

if len(sys.argv) == 3:
  output = sys.argv[2]
else:
  output = os.path.splitext(filename)[0]+'.pbin'

A = scipy.io.mmread(filename)
io = PBIO.PetscBinaryIO()
io.writeMatSciPy(open(output,'w'),A) 






