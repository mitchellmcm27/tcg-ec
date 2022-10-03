#!/usr/bin/env bash

# quicky shell script to build the fo_h2O system

# clean out any previous builds
echo cleaning prevous builds
rm database/*.tar.gz
rm endmembers/*.emml
rm phases/*.phml
rm reactions/*.rxml
rm -rf reactions/fo_h2O_hydration

# generate spudfiles and build  thermodynamic database
jupyter nbconvert --to notebook --execute notebooks/Generate_berman_endmembers.ipynb
jupyter nbconvert --to notebook --execute notebooks/Generate_phases.ipynb
tcg_builddb --just_src --include_swim -zi database --calibfile calib_params.csv

#generate spudfiles and libraries for reactions
jupyter nbconvert --to notebook --execute notebooks/Generate_reactions.ipynb
cd reactions
tcg_buildrx *.rxml --include_swim -i
module load ./fo_h2o_hydration/fo_h2o_hydration.module

# test
cd ../tests
pytest  --disable-warnings test*
cd ..

# clean up
rm -f notebooks/*.nbconvert.ipynb
rm -rf *.build reactions/*.build
module unload reactions/fo_h2o_hydration/*.module

exit 0



