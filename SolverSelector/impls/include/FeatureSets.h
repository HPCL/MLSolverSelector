// These are all the feature sets that can be used. To define a new feature set, copy 
// the full feature set ( in the else statement, put it in a elif statement, and then 
// delete any features that are not needed. Then in the makefile create a new build
// for the new feature set name. Basically, you need to set FSET through the makefile.   
#ifdef FULL_DEBUGFAIL
  
  // Sum based features 
  #define NNZ _NNZ 
  #define AVGNONZEROSPERROW _AVGNONZEROSPERROW 
  #define ABSOLUTENONZEROSUM _ABSOLUTENONZEROSUM 
  #define FROBENIUSNORM _FROBENIUSNORM 
  #define AVERAGEDIAGONALDISTANCE _AVERAGEDIAGONALDISTANCE            
  #define ROWVARIANCE _ROWVARIANCE 
  #define TRACE _TRACE 
  #define ABSOLUTETRACE _ABSOLUTETRACE 
  #define DIAGONALMEAN _DIAGONALMEAN 
  #define DIAGONALAVERAGE _DIAGONALAVERAGE    
  #define DIAGONALNONZEROS _DIAGONALNONZEROS   
  #define SUMCOUNT (NNZ + AVGNONZEROSPERROW  + ABSOLUTENONZEROSUM + FROBENIUSNORM  + \
      AVERAGEDIAGONALDISTANCE + ROWVARIANCE + TRACE + ABSOLUTETRACE + DIAGONALNONZEROS +DIAGONALMEAN +  DIAGONALAVERAGE ) 

  // max based features 
  #define INFINITYNORM _INFINITYNORM 
  #define MAXNONZEROSPERROW _MAXNONZEROSPERROW 
  #define DIAGONALSIGN _DIAGONALSIGN 
  #define LOWERBANDWIDTH _LOWERBANDWIDTH 
  #define UPPERBANDWIDTH _UPPERBANDWIDTH 
  #define ROWDIAGONALDOMINANCE _ROWDIAGONALDOMINANCE 
  #define MAXCOUNT (INFINITYNORM + MAXNONZEROSPERROW + DIAGONALSIGN + LOWERBANDWIDTH + UPPERBANDWIDTH + ROWDIAGONALDOMINANCE )
  
  //min based features
  #define MINNONZEROSPERROW _MINNONZEROSPERROW 
  #define MINCOUNT (MINNONZEROSPERROW)
  
  //column based features ( need a column of values to be communicated for calculation on root) 
  #define ONENORM _ONENORM 
  #define COLDIAGONALDOMINANCE _COLDIAGONALDOMINANCE 
  #define COLUMNVARIANCE _COLUMNVARIANCE 
  #define COLCOUNT (ONENORM + COLDIAGONALDOMINANCE + COLUMNVARIANCE )
  

  // symmetry based features ( need all sample data to be communicated for calculation on root) 
  #define SYMMETRICITY _SYMMETRICITY 
  #define NONZEROPATTERNSYMMETRY _NONZEROPATTERNSYMMETRY 
  #define SYMMETRICINFINITYNORM _SYMMETRICINFINITYNORM 
  #define SYMMETRICFROBENIUSNORM _SYMMETRICFROBENIUSNORM 
  #define ANTISYMMETRICINFINITYNORM _ANTISYMMETRICINFINITYNORM 
  #define ANTISYMMETRICFROBENIUSNORM _ANTISYMMETRICFROBENIUSNORM
  #define SYMMETRY 1 

#endif


#ifdef RS1_DEBUGFAIL
  
  // Sum based features 
  #define NNZ _NNZ
  #define ABSOLUTETRACE _ABSOLUTETRACE 
  #define DIAGONALAVERAGE _DIAGONALAVERAGE    
  #define DIAGONALNONZEROS _DIAGONALNONZEROS   
  #define SUMCOUNT ( NNZ + ABSOLUTETRACE + DIAGONALNONZEROS +  DIAGONALAVERAGE ) 
  
  // Max based features
  #define INFINITYNORM _INFINITYNORM
  #define LOWERBANDWIDTH _LOWERBANDWIDTH 
  #define MAXCOUNT (INFINITYNORM + LOWERBANDWIDTH )
  
  // Min based features 
  #define MINNONZEROSPERROW _MINNONZEROSPERROW
  #define MINCOUNT (MINNONZEROSPERROW)
  
  // Column based features 
  #define COLUMNVARIANCE _COLUMNVARIANCE 
  #define COLCOUNT COLUMNVARIANCE
  
  // Symmetry based features  
  #define NONZEROPATTERNSYMMETRY _NONZEROPATTERNSYMMETRY
  #define SYMMETRY 1 
  
#endif 

#ifdef RS2_DEBUGFAIL

  #define NNZ _NNZ
  #define ABSOLUTETRACE _ABSOLUTETRACE 
  #define DIAGONALAVERAGE _DIAGONALAVERAGE    
  #define DIAGONALNONZEROS  _DIAGONALNONZEROS   
  #define SUMCOUNT (NNZ + ABSOLUTETRACE + DIAGONALNONZEROS + DIAGONALAVERAGE ) 

  #define LOWERBANDWIDTH _LOWERBANDWIDTH 
  #define MAXCOUNT LOWERBANDWIDTH 
  
  #define MINNONZEROSPERROW _MINNONZEROSPERROW
  #define MINCOUNT MINNONZEROSPERROW
  
  #define COLUMNVARIANCE _COLUMNVARIANCE
  #define COLCOUNT COLUMNVARIANCE 

  #define SYMMETRY 0 


#endif

//DEBUGGING TODO FIXME AAAAAH 
#define RS1_KANIKA 

#ifdef RS1_KANIKA

  #define ABSOLUTETRACE _ABSOLUTETRACE
  #define TRACE _TRACE 
  #define DIAGONALMEAN _DIAGONALMEAN 
  #define DIAGONALNONZEROS _DIAGONALNONZEROS   
  #define SUMCOUNT (TRACE + ABSOLUTETRACE + DIAGONALNONZEROS +DIAGONALMEAN ) 


  #define ONENORM _ONENORM 
  #define COLCOUNT (ONENORM)
  
  #define MAXCOUNT 0
  #define MINCOUNT 0
  
  #define SYMMETRICINFINITYNORM _SYMMETRICINFINITYNORM 
  #define SYMMETRY 1 
#endif


#ifdef RSKANIKA 
  
  #define AVGNONZEROSPERROW _AVGNONZEROSPERROW  
  #define AVERAGEDIAGONALDISTANCE _AVERAGEDIAGONALDISTANCE                 
  #define TRACE _TRACE
  #define DIAGONALAVERAGE _DIAGONALAVERAGE     
  #define DIAGONALNONZEROS  _DIAGONALNONZEROS 
  #define SUMCOUNT ( AVGNONZEROSPERROW +AVERAGEDIAGONALDISTANCE + TRACE + DIAGONALAVERAGE + DIAGONALNONZEROS ) 


  #define MAXCOUNT 0
  #define MINCOUNT 0

  #define COLDIAGONALDOMINANCE _COLDIAGONALDOMINANCE 
  #define COLCOUNT COLDIAGONALDOMINANCE 
  
  #define SYMMETRICFROBENIUSNORM _SYMMETRICFROBENIUSNORM
  #define SYMMETRY 1 

#endif

