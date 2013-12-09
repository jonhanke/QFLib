#!/bin/bash

####################################
##  Install script for Pari-2.1.7 ##
####################################

## Uncompress the Tarball
tar -zxvf ../Tarballs/pari-2.1.7.tgz -C ../Builds/

## Build the software
cd ../Builds/pari-2.1.7
./Configure  --prefix=../..
make all
make bench

## Install the software
make install