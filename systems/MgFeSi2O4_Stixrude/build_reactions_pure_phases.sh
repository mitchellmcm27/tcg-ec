#!/usr/bin/env bash

# quicky shell script to build the Mg2SIO4 Stixrude  system

# clean out any previous builds
echo cleaning prevous builds
ROOT=Mg2SiO4_stixrude
rm reactions/${ROOT}.rxml
rm -rf reactions/${ROOT}

if [ ! -f database/MgFeSi2O4_Stixrude.tar.gz ]; then
    echo 'rebuilding database'
    ./build_database.sh
fi


#generate spudfiles and libraries for reactions
jupyter nbconvert --to notebook --execute notebooks/Generate_reactions_Mg_pure_phases.ipynb
cd reactions
tcg_buildrx ${ROOT}.rxml -i
module load ./${ROOT}/${ROOT}.module

# test
cd ../tests
pytest  --disable-warnings test*_pure*
cd ..


# clean up
rm -f notebooks/*.nbconvert.ipynb
rm -rf *.build reactions/*.build

exit 0