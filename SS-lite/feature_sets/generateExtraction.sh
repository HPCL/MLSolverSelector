
#!/usr/bin/env bash

#### This script generates the source code required for a new featureExtraction algorithm. Using this, you
### can generate a new class for each new feature set. 

### All you need to do is provide the name of the feature set file containing the appropraite defines (see 
### full_feature_set.fs and example1.fs and the name of the class you want to create. 

### Once you generate a class, the script will insert a if statement into the appropriate place in the 
### source code so that you can then use your feature set for feature extraction with the command line 
### arguement --feature_set <NAME>  -- After doing this you must recompile. 

#### Step One -- Cat the files together. 

### Call with ./generateExtraction.sh <NAME> <featureExtraction_file> 

## First, clear the source directory 
rm -f include/*.h

cp templates/PetscFeatureSetSelector.template include/PetscFeatureSetSelector.h


for i in impls/*.fs; do
    [ -f "$i" ] || break
    fname=`basename -s .fs "$i"`
    cat templates/FeatureExtractionHead.template $i templates/FeatureExtractionTail.template | sed -e 's/__NAME__/'"$fname"'/g' > include/${fname}_featureExtraction.h
    sed -i -e 's@__SEARCH__IF@Feature Set \n    if ( featureSet == \"'"$fname"'\" ) return FeatureExtraction_'"$fname"'::ExtractJacobianFeatures(J,edge,interior,fnames,matvecs); \n\n    \/\/__SEARCH__IF@g' include/PetscFeatureSetSelector.h 
    sed -i -e 's@__SEARCH__INCLUDE@ Feature Set Include \n#include \"'"$fname"'_featureExtraction.h\" \n\n\/\/__SEARCH__INCLUDE@g' include/PetscFeatureSetSelector.h 
    sed -i -e 's@__SEARCH__EXISTS@   Feature Set Exists \n   if( featureSet == \"'"$fname"'\") return true; \n\n\/\/__SEARCH__EXISTS@g' include/PetscFeatureSetSelector.h 
done


#cat FeatureExtractHead.h $2 FeatureExtractTail.h | sed -e 's/__NAME__/'"$1"'/g' > ${1}_featureExtraction.h

#sed -i -e 's@__SEARCH__AND__REPLACE__TAG@Feature Set \n    if ( featureSet == \"'"$1"'\" ) return FeatureExtraction'"$1"'_name::ExtractJacobianFeatures(J,edge,interior,fnames,matvecs); \n\n    \/\/__SEARCH__AND__REPLACE__TAG@g' PetscFeatureSetSelector.hpp 
