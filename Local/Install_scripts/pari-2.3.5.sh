#!/bin/bash

####################################
##  Install script for Pari-2.3.5 ##
####################################

## Uncompress the Tarball
tar -zxvf ../Tarballs/pari-2.3.5.tgz -C ../Builds/

## Build the software
cd ../Builds/pari-2.3.5
./Configure  --prefix=../..
make all
make bench

## Install the software
make install