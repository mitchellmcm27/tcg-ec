#!/usr/bin/env bash

# quicky shell script to build phases in the MgFeSi2O4 Stixrude  system

# clean out any previous builds
echo cleaning prevous builds
rm database/*.tar.gz
rm endmembers/*.emml
rm phases/*.phml

# generate spudfiles and build  thermodynamic database
jupyter nbconvert --to notebook --execute notebooks/Generate_stixrude_endmembers.ipynb
jupyter nbconvert --to notebook --execute notebooks/Generate_pure_phases.ipynb
jupyter nbconvert --to notebook --execute notebooks/Generate_solution_phases.ipynb

tcg_builddb --just_src  -zi database

exit 0



