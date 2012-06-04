#!/bin/bash
#
# Set the name of the job.
#$ -N main_1000-1010
#
# Make sure that the .e and .o file arrive in the 
# working directory
#$ -cwd
# Merge the standard out and standard error to one file
#$ -j y 
#

## This is my job! =)
./main 1000 1010
