#!/bin/bash

#########################################################
## Script to erase all files in the "Local" directory  ##
## (except Install_scripts and Tarballs),              ##
#########################################################

## Remove all old files (except tarballs and install scripts)
rm -rf ../Builds/*
rm -rf ../include/*
rm -rf ../docs/*
rm -rf ../bin/*
rm -rf ../lib/*
rm -rf ../man/*
rm -rf ../share/*
