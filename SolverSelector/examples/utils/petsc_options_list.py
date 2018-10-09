##### This is a script that gets the solver information, preconditioner information, and all the
##### parameters from petsc

import re
import glob
import os
from subprocess import call 

###TODO Add inputs for solvers and preconditioners, and hypre preconditioners 
### to make a config file only containing certain solvers 

##### Step one is two call the scanner function 
call(["make", "petsc_scan.a"])

#### Step 2, get a list of all the solvers and preconditioners
f = open("petscoutput1","w")
call(["./petsc_scan.a", "-pc_type", "hypre", "-help"], stdout=f)
f.close()

f = open("petscoutput1","r")
solvers = []
precond = []
hypre = []
ksp_pattern = re.compile("\-ksp_type \<.*?\>: Krylov method \(one of\) (.*?) \(KSPSetType\)")
pc_pattern = re.compile( "\-pc_type \<.*?\>: Preconditioner \(one of\) (.*?) \(PCSetType\)")
hy_pattern = re.compile( "\-pc_hypre_type \<.*?\> \(choose one of\) (.*?) \(")
for line in f:
   match_ksp = ksp_pattern.findall(line)
   if ( match_ksp ):
        solvers = match_ksp[0].split(" ")

   match_pc = pc_pattern.findall(line)
   if match_pc :
            precond = match_pc[0].split(" ")
   
   match_hy = hy_pattern.findall(line)
   if match_hy:
            hypre = match_hy[0].split(" ")
            
#### Step 3, get a list of standard ksp solver options 
petsc_ignore = ["monitor","plot"] ###useless options for SS 


gsolver_dict = {}
for ksp in solvers:
    gsolver_dict[ksp] = { "Parameters" : [] , "Description" : [] }
    f = open("petscoutput1","w")
    call(["./petsc_scan.a","-ksp_type",ksp,"-help"], stdout=f)
    f.close()
    
    
    ### step one, get all the solver parameters ( could be nested )
    for ksp1 in solvers:
        ss = "(\-ksp\_"+ksp1+"\_[a-z\_A-Z]*)(.*)"
        p1 = re.compile(ss)
    
        f = open("petscoutput1","r")         
        for line in f:
            match = p1.findall(line)
            if (match):
                gsolver_dict[ksp]["Parameters"].append(match[0][0])
                gsolver_dict[ksp]["Description"].append(match[0][1])
        f.close()        

gprecond_dict = {}
for pc in precond:
    if pc != "hypre":
        
        gprecond_dict[pc] =  { "Parameters" : [] , "Description" : [] }
        f = open("petscoutput1","w")
        call(["./petsc_scan.a","-pc_type",pc,"-help"], stdout=f)
        f.close()
        
        ss = "(\-pc\_[a-z\_A-Z]*)(.*)"
        p1 = re.compile(ss)
        
        f = open("petscoutput1","r")         
        for line in f:
            match = p1.findall(line)
            if (match):
                if  match[0][0] != "-pc_type":
                    gprecond_dict[pc]["Parameters"].append(match[0][0])
                    gprecond_dict[pc]["Description"].append(match[0][1])

        f.close()
        
### Get the hypre info
ghypre_dict = {}
for htype in hypre:
    ghypre_dict[htype] =  { "Parameters" : [] , "Description" : [] }
    f = open("petscoutput1","w")
    call(["./petsc_scan.a","-pc_type","hypre","-pc_hypre_type",htype,"-help"], stdout=f)
    f.close()
        
    ss = "(\-pc\_[a-z\_A-Z]*)(.*)"
    p1 = re.compile(ss)
        
    f = open("petscoutput1","r")         
    for line in f:
        match = p1.findall(line)
        if (match):
           if  match[0][0] != "-pc_type" and match[0][0] != "-pc_hypre_type" :
              ghypre_dict[htype]["Parameters"].append(match[0][0])
              ghypre_dict[htype]["Description"].append(match[0][1])

    f.close()


def print_solver( file_name, basename, parameters, description, types ):
  if basename != "python": 

    def_pat = re.compile("\<(.*?)\>")  
    if ( types == "solver" ):
        file_name.write("@SOLVER %s \n" %basename )
    elif ( types == "precond" ):            
        file_name.write("@PRECON %s \n" %basename)
    elif ( types == "hypre"):            
        file_name.write("@PRECON hypre %s \n" % basename )

    for i,p in enumerate(parameters):                   
         
         def_value = def_pat.findall(description[i]) 
         def1 = "none"
         if def_value:
            def1 = def_value[0]   
         
         file_name.write( "@PARAMETER %s %s   # %s \n" %(p,def1, description[i])  );  
  file_name.write("\n")
  return 


file_name = open("PetscSolverOptions.db","w")

file_name.write( "# List of options for each solver \n\n" );

for solver in sorted(gsolver_dict.keys()):     
    print_solver( file_name, solver, gsolver_dict[solver]["Parameters"], gsolver_dict[solver]["Description"], "solver" )
file_name.close()

file_name = open("PetscPreconOptions.db","w")
file_name.write( "# List of options for each preconditioner \n\n" );

for pc in sorted(gprecond_dict.keys()):
    print_solver( file_name, pc, gprecond_dict[pc]["Parameters"], gprecond_dict[pc]["Description"], "precond" );

for hy in sorted(ghypre_dict.keys()):
    print_solver( file_name, hy, ghypre_dict[hy]["Parameters"], ghypre_dict[hy]["Description"], "hypre" );


file_name.close()
