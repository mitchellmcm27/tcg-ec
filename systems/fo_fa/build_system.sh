#!/usr/bin/env bash

# quicky shell script to build the fo_fa  system

# clean out any previous builds
echo cleaning prevous builds
rm database/*.tar.gz
rm endmembers/*.emml
rm phases/*.phml
rm reactions/fo_fa_binary.rxml
rm -rf reactions/fo_fa_binary

# generate spudfiles and build  thermodynamic database
jupyter nbconvert --to notebook --execute notebooks/Generate_berman_endmembers.ipynb
jupyter nbconvert --to notebook --execute notebooks/Generate_xmelts_endmembers.ipynb
jupyter nbconvert --to notebook --execute notebooks/Generate_phases.ipynb
tcg_builddb --just_src  -zi database

#generate spudfiles and libraries for reactions
jupyter nbconvert --to notebook --execute notebooks/Generate_reactions.ipynb
cd reactions
tcg_buildrx fo_fa_binary.rxml -i
module load ./fo_fa_binary/fo_fa_binary.module

# test
cd ../tests
pytest  --disable-warnings test*
cd ..


# clean up
rm -f notebooks/*.nbconvert.ipynb
rm -rf *.build reactions/*.build
module unload reactions/fo_fa_binary/*.module

exit 0



