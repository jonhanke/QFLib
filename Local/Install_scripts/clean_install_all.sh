#!/bin/bash

#########################################################
## Script to erase all files in the "Local" directory  ##
## (except Install_scripts and Tarballs),              ##
## and run all install scripts                         ##
#########################################################

## Remove all old files (except tarballs and install scripts)
./delete_all.sh

## Create a Local/Builds directory to hold the dependency builds
mkdir ../Builds

## Rebuild GMP
./gmp-5.0.2.sh

## Rebuild Pari
#./pari-2.5.0.sh     ## This builds, but is missing some of the header declarations I used in my older C++ code.
#./pari-2.2.10.alpha.sh    ## This fails to build, and is sensitive to spaces in the global path! =(
#./pari-2.1.7.sh      ## This has some readline error with another installed version in /usr/local. =(
./pari-2.3.5.sh
