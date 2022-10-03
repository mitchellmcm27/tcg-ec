#!/usr/bin/env bash

# quicky shell script to build a clean fo_si02 database and reactions and run tests
export REACTIONS=fo_sio2_poly_linear_rxns

# clean out any previous builds
echo cleaning prevous builds
rm -f database/*.tar.gz
rm -f endmembers/*.emml
rm -f phases/*.phml
rm -f reactions/*.rxml
rm -rf reactions/${REACTIONS}

# generate spudfiles and build  thermodynamic database
cd notebooks
for d in endmembers phases
do
    files=$d/*.ipynb
    for file in $files
    do
        jupyter nbconvert --to notebook --execute $file
    done
done
cd ..
tcg_builddb --just_src  -l fo_sio2_db -zi database


#generate spudfiles and libraries for reactions
jupyter nbconvert --to notebook --execute notebooks/reactions/Generate_${REACTIONS}.ipynb
cd reactions
tcg_buildrx *.rxml  -i
module load ./${REACTIONS}/*.module
module list
cd ..

# clean up
rm -rf *.build reactions/*.build
rm -f notebooks/*/*.nbconvert.ipynb

# test
cd tests
pytest  --disable-warnings test-*

cd ../reactions
module unload ${REACTIONS}/*.module
exit 0





