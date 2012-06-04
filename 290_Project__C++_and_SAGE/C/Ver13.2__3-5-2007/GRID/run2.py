#!/usr/bin/python

import os
import time



# This is a function that takes a list of form numbers, and computes those forms! =)
def compute_aux_quaternaries(form_number_list):
    """Takes a list of form numbers, and computes the numbers represented by them"""

    Project_Name = "290_Project__Auxiliary_Quaternary_Escalators"
    Small_Prime_File = "primelist_upto_15millionth_prime.txt"
    Big_Prime_File = "primelist_upto_40millionth_prime.txt"
    Form_File = "104_auxiliary_quaternaries.txt"
    Cusp_Dir = "Auxiliary_Constants"
    Cusp_Prefix = "Aux_const_"

    Output_Log_Dir = "Aux_C++_Output_Logs"
    Command_start = "main " + Project_Name + " "
    Command_middle = " " + Form_File + " " + Cusp_Dir + " " + Cusp_Prefix + " "
    
    
    for form_num in form_number_list:

        print
        print " About to compute form #", form_num 
        print time.ctime()

        command_line = Command_start + Small_Prime_File + Command_middle + str(form_num) + " > " + Output_Log_Dir + "/Form" + str(form_num) + ".txt"
        print command_line
        os.system(command_line)
        
#        if (form_num != 6414):
#            os.system(Command_start + Small_Prime_File + Command_middle + str(form_num) + " > " + Output_Log_Dir + "/Form" + str(form_num) + ".txt")
#        else:
#            os.system(Command_start + Big_Prime_File + Command_middle + str(form_num) + " > " + Output_Log_Dir + "/Form" + str(form_num) + ".txt")

            
    print time.ctime()



# This looks through each of the Form output files to see if there is an error
def find_form_errors():

    filename_base_string = "OUTPUT_LOGS/Form"

    Form_Error_list = []

    for num in range(1,6561):

        # Check if the output file exists
        form_filename = filename_base_string + str(num) + ".txt"
        if os.path.exists(form_filename):
            
            # Read in the form output file
            form_file = open(form_filename, "r")
            output_lines = form_file.readlines()
            form_file.close()

            # Set the error flag
            error_flag = false
                
            # Check for explicit errors
            for line in output_lines:
                if (line.count("Error") > 0) and (line.count("_ReadQuaternaryTheta() Error:  The precision of the theta function is not 10,000... =( ... so we'll have to regenerate it...") == 0):
                    error_flag = true


            # Check that the local cover worked


            # Check for a final report of square-free exceptions
            

            # Check for a program ending

            
                    

            # Report an error
            print "\nError for Form", num       
            Form_Error_list.append(num)


    # Return the list of forms with an error!
    return Form_Error_list



# This deletes a range of jobs (by Job ID)
def qdel(JobId_list): 
    for num in JobId_list: 
        os.system("qdel " + str(num))    


###################################################################################################
# This is a function that takes a list of form numbers, and computes those forms on the grid! =) ##
###################################################################################################
def compute_aux_quaternaries__on_grid(form_number_list):
    """Takes a list of form numbers, and computes the numbers represented by them using the grid"""

    Project_Name = "290_Project__Auxiliary_Quaternary_Escalators"
    Small_Prime_File = "primelist_upto_15millionth_prime.txt"
    Big_Prime_File = "primelist_upto_40millionth_prime.txt"
    Form_File = "104_auxiliary_quaternaries.txt"
    Cusp_Dir = "Auxiliary_Constants"
    Cusp_Prefix = "Aux_const_"
    Output_Log_Dir = "Aux_C++_Output_Logs"
    
    Job_Name = "_290_Project__104_Auxiliaries___Form"

      
    for form_num in form_number_list:

        # Make the qsub file for this form
        print "Make the SGE job file for sending the Form #"+ str(form_num) + " computation to the grid"
    
        jobfile_template = """
#!/bin/bash
#
## Use the BASH shell for this script
#$ -S /bin/bash


## Set the name of the job  (i.e. $JOB_NAME)
#$ -N NNNNN


## Set the output file path
#$ -o JJJJJ


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


export LD_LIBRARY_PATH=/usr/local/lib:/home/postdoc/jonhanke/bin/Installed_Software/LOCAL/lib

echo " Submitting job #$SGE_TASK_ID at $(date)."
echo
echo ' This is task #'$SGE_TASK_ID
./main PROJECT_NAME PRIME_FILE FORM_FILE CUSP_DIR CUSP_PREFIX XXXX >> OUTPUT_LOG_DIR/FormXXXX.txt
echo
echo " Finished submitting job #$SGE_TASK_ID at $(date)."
"""

        # Customize the jobfile text for the current form
        jobfile_text = jobfile_template
        jobname = "_290_basic_quaternary_output__Form" + str(form_num)
        jobfile_text = jobfile_text.replace('NNNNN', Job_Name)
        jobfile_text = jobfile_text.replace('JJJJJ', jobname)
        jobfile_text = jobfile_text.replace('XXXX', str(form_num))
        jobfile_text = jobfile_text.replace('PROJECT_NAME', Project_Name)
        jobfile_text = jobfile_text.replace('PRIME_FILE', Small_Prime_File)
        jobfile_text = jobfile_text.replace('FORM_FILE', Form_File)
        jobfile_text = jobfile_text.replace('CUSP_DIR', Cusp_Dir)
        jobfile_text = jobfile_text.replace('CUSP_PREFIX', Cusp_Prefix)
        jobfile_text = jobfile_text.replace('OUTPUT_LOG_DIR', Output_Log_Dir)
        
        # Writing the SGE job file
        jobfile_name = "jobfile_for_form" + str(form_num)
        jobfile = open(jobfile_name, "w")
        jobfile.write(jobfile_text)
        jobfile.close()


        # Run qsub on the jobfile
        print "Queueing the job for Form #" + str(form_num)
        os.system("qsub " + jobfile_name)

        # Remove the jobfile
        time.sleep(0.1)
#        os.system("rm -f " + jobfile_name)




###################################################################################################
# This is a function that takes a list of form numbers, and computes those forms on the grid! =) ##
###################################################################################################
def compute_basic_quaternaries__on_grid(form_number_list):
    """Takes a list of form numbers, and computes the numbers represented by them using the grid"""


    #Project_Name = "290_Project__Basic_Quaternary_Escalators__8-21-05__with_Approx_Thetas"
    Project_Name = "290_Project__Basic_Quaternary_Escalators__11-2-06"
    Small_Prime_File = "primelist_upto_15millionth_prime.txt";
    Big_Prime_File = "primelist_upto_40millionth_prime.txt";
    Job_Name = "_290_Project__6560_Quaternaries___Form"
    Form_File = "290_Project__Basic_Quaternary_Escalators__form_list.txt"
    Cusp_Dir = "Basic_Constants"
    Cusp_Prefix = "Basic_const_"
    Output_Log_Dir = "Basic_C++_Output_Logs"
      


    for form_num in form_number_list:

        # Make the qsub file for this form
        print "Make the SGE job file for sending the Form #"+ str(form_num) + " computation to the grid"
    
        jobfile_template = """
#!/bin/bash
#
## Use the BASH shell for this script
#$ -S /bin/bash


## Set the name of the job  (i.e. $JOB_NAME)
#$ -N NNNNN


## Set the output file path
#$ -o JJJJJ


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


export LD_LIBRARY_PATH=/usr/local/lib:/home/postdoc/jonhanke/bin/Installed_Software/LOCAL/lib

echo " Submitting job #$SGE_TASK_ID at $(date)."
echo
echo ' This is task #'$SGE_TASK_ID
../main PROJECT_NAME PRIME_FILE FORM_FILE CUSP_DIR CUSP_PREFIX XXXX >> OUTPUT_LOG_DIR/FormXXXX.txt
echo
echo " Finished submitting job #$SGE_TASK_ID at $(date)."
"""

        # Customize the jobfile text for the current form
        jobfile_text = jobfile_template
        jobname = "_290_basic_quaternary_output__Form" + str(form_num)
        jobfile_text = jobfile_text.replace('NNNNN', Job_Name)
        jobfile_text = jobfile_text.replace('JJJJJ', jobname)
        jobfile_text = jobfile_text.replace('XXXX', str(form_num))
        jobfile_text = jobfile_text.replace('PROJECT_NAME', Project_Name)
        if (form_num != 6414):
            jobfile_text = jobfile_text.replace('PRIME_FILE', Small_Prime_File)
        else:
            jobfile_text = jobfile_text.replace('PRIME_FILE', Big_Prime_File)
        jobfile_text = jobfile_text.replace('FORM_FILE', Form_File)
        jobfile_text = jobfile_text.replace('CUSP_DIR', Cusp_Dir)
        jobfile_text = jobfile_text.replace('CUSP_PREFIX', Cusp_Prefix)
        jobfile_text = jobfile_text.replace('OUTPUT_LOG_DIR', Output_Log_Dir)
        
        
        # Writing the SGE job file
        jobfile_name = "jobfile_for_form" + str(form_num)
        jobfile = open(jobfile_name, "w")
        jobfile.write(jobfile_text)
        jobfile.close()


        # Run qsub on the jobfile
        print "Queueing the job for Form #" + str(form_num)
        os.system("qsub " + jobfile_name)

        # Remove the jobfile
        time.sleep(0.1)
#        os.system("rm -f " + jobfile_name)

