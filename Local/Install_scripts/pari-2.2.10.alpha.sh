#!/bin/bash

###########################################
##  Install script for Pari-2.2.10.alpha ##
###########################################

## Uncompress the Tarball
tar -zxvf ../Tarballs/pari-2.2.10.alpha.tgz -C ../Builds/

## Build the software
cd ../Builds/pari-2.2.10.alpha
./Configure  --prefix=../..
make all
make bench

## Install the software
make install