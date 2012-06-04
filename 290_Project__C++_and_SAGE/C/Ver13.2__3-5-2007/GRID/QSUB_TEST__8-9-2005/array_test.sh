#!/bin/bash
#
## Use the BASH shell for this script
#$ -S /bin/bash


## Set the name of the job  (i.e. $JOB_NAME)
#$ -N _290_Project__6560_Quaternaries___Form



##################################################################################
### THIS IS TOTALLY USELESS -- Can't get a bash variable into a qsub variable! ###
##################################################################################
#SGE_JOB_NAME="$(date +%F__%H:%M)__12"
#JOB_NAME="$(date +%F__%H:%M)__23"
#export SGE_JOB_NAME="$(date +%F__%H:%M)"
#echo $SGE_JOB_NAME
#set TEST_OUT = "XXX_290_Project--TEST__$(date +%F__%H:%M)__Form\#$SGE_TASK_ID"
##################################################################################



## Set the output file path
#$ -o $JOB_NAME_$TASK_ID


## Send me mail when a job aborts
#$ -M jonhanke@math.duke.edu
### -m e


# Make sure that the .e and .o file arrive in the 
#working directory
#$ -cwd
#
# Merge the standard out and standard error to one file
#$ -j y 
#
# My code is re-runnable
#$ -r y


echo " Submitting job #$SGE_TASK_ID at $(date)."
echo
echo ' This tests the array feature of qsub!'
echo ' This is task #'$SGE_TASK_ID
echo
echo " Finished submitting job #$SGE_TASK_ID at $(date)."
