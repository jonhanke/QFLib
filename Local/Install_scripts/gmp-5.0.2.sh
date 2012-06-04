#!/bin/bash

####################################
##  Install script for Pari-2.1.7 ##
####################################

## Uncompress the Tarball
tar -zxvf ../Tarballs/gmp-5.0.2.tar.bz2 -C ../Builds/

## Get the absolute (dereferenced) pathname of the "Local" directory
cd ..
Local_Dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## Build the software
cd Builds/gmp-5.0.2
./configure  --prefix=$Local_Dir  --enable-cxx
make 
make check

## Install the software
make install