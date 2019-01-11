 
  // HERE WE NEED TO SPECIFY THE FEATURE SET 
  // First step is to copy the full feature set from the feature_set directory. 
  //
  // The names files suggests the fset is AbsoluteTrace, OneNorm, ColDiagDom, SymmetricInf, DiagonalMean, Frob, nnz 
  //
  // So, we leave those one and remove all others. Anything you dont want, just set 0 and comment out the rest 


   // To install a feature set you need to copy the .fs file into the feature_sets/impls directory 
   // and run the ./generate... script -- then recompile. 

 
  #define ABSOLUTETRACE _ABSOLUTETRACE 
  #define ONENORM _ONENORM 
  #define COLDIAGONALDOMINANCE _COLDIAGONALDOMINANCE 
  #define SYMMETRICINFINITYNORM _SYMMETRICINFINITYNORM 
  #define DIAGONALMEAN _DIAGONALMEAN 
  #define FROBENIUSNORM _FROBENIUSNORM 
  #define NNZ _NNZ 

  // Everything else gets defined to zero ( in Vim :%s/_/ 0 \/\/ _/gc  then n,n,n, ...., y,y,y,y,y,y )  

  #define AVGNONZEROSPERROW 0 // _AVGNONZEROSPERROW 
  #define ABSOLUTENONZEROSUM 0 //  _ABSOLUTENONZEROSUM 
  #define AVERAGEDIAGONALDISTANCE  0 // _AVERAGEDIAGONALDISTANCE            
  #define ROWVARIANCE  0 // _ROWVARIANCE 
  #define TRACE  0 // _TRACE 
  #define DIAGONALAVERAGE  0 // _DIAGONALAVERAGE    
  #define DIAGONALNONZEROS  0 // _DIAGONALNONZEROS   
  #define INFINITYNORM  0 // _INFINITYNORM 
  #define MAXNONZEROSPERROW  0 // _MAXNONZEROSPERROW 
  #define DIAGONALSIGN  0 // _DIAGONALSIGN 
  #define LOWERBANDWIDTH  0 // _LOWERBANDWIDTH 
  #define UPPERBANDWIDTH  0 // _UPPERBANDWIDTH 
  #define ROWDIAGONALDOMINANCE  0 // _ROWDIAGONALDOMINANCE 
  #define MINNONZEROSPERROW  0 // _MINNONZEROSPERROW 
  #define COLUMNVARIANCE  0 // _COLUMNVARIANCE 
  #define SYMMETRICITY  0 // _SYMMETRICITY 
  #define NONZEROPATTERNSYMMETRY  0 // _NONZEROPATTERNSYMMETRY 
  #define SYMMETRICFROBENIUSNORM  0 // _SYMMETRICFROBENIUSNORM 
  #define ANTISYMMETRICFROBENIUSNORM  0 // _ANTISYMMETRICFROBENIUSNORM
  #define ANTISYMMETRICINFINITYNORM  0 // _ANTISYMMETRICINFINITYNORM 


