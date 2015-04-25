# QFLIB Installation Notes #

The following steps should allow you to install QFLIB on Linux:

  1. Clone the repository from Google code by typing
```
    hg clone https://code.google.com/p/qflib/
```
> > in the desired directory.  This will make a `qflib` subdirectory with the code repository.
  1. Install the dependencies (GMP and Pari) by typing
```
    cd qflib/Local/Install_scripts
    ./clean_install_all.sh
```
> > This will build the required software in `qflib/Local/Builds` and install it in `qflib/Local`.
  1. Set the appropriate paths in the `Makefile`, via
```
    cd ../../qflib/290_Project__C++_and_SAGE/C/Current_Version
    emacs Makefile
```
  1. To build the `main.cc` C code and associated Doxygen documentation, type
```
    make
```
  1. Set some required environment variables
```
    export QFLIB_local_project_data_dir='/tmp/QFLIB_Project_Data/'
    export QFLIB_remote_user_at_machine='me@some_machine.com'
    export QFLIB_remote_MAGMA_path='/usr/local/bin/magma'
    export QFLIB_remote_temp_dir='/tmp/QFLIB_temporary_computations/'
```
  1. Set an optional environment variables if we want to relocate the tables of primes files
```
    export QFLIB_local_primes_dir='../../Primes/'      // This is the default setting
```
  1. Setup an ssh key-pair between the machine running the QFLIB software and the remote account performing the cuspidal computations (specified by the shell variable `QFLIB_remote_user_at_machine`).  If the local machine is set for these computations (e.g. `user@localhost`) then this step is unnecessary.


---


# Misc QFLIB Notes #

  1. All QFLIB project files are stored in a subfolder of the folder specified by the shell variable `QFLIB_local_project_data_dir`.
  1. Sample code to check representability
```
./main test_project blah 104_auxiliary_quaternaries.txt  "" cusp_const__ 10
```


# Warnings/To Do #

  1. There are still many hardcoded paths... e.g. in the Makefile.  These need to be converted to use shell variable settings!


---


# 290-Theorem Notes #
  1. To download the cupsidal constant datafiles used for the 290-Theorem use
```
    mkdir ~/290-Theorem_cusp_constants
    cd ~/290-Theorem_cusp_constants
    wget http://data.jonhanke.com/290_Data_Tarballs/Basic_and_Auxiliary_Constants.tgz
    tar -zxvf Basic_and_Auxiliary_Constants.tgz
```
> > This will download the datafiles into two directories `Basic_Constants` and `Auxiliary_Constants` whose filenames are of the form `Basic_const_1004.txt` or `Aux_const_45.txt`.  Therefore the associated cusp\_const\_prefix to use for these files are `Basic_const_` and `Aux_const_` respectively.
  1. To run the 290-Theorem computation on these files we change back to the QFLIB main executable directory and use the commands
```
    ./main basic_project blah 6560_basic_quaternaries.txt ~/290-Theorem_cusp_constants/Basic_Constants/ Basic_const_ 1 6560
    ./main aux_project blah 104_auxiliary_quaternaries.txt  ~/290-Theorem_cusp_constants/Auxiliary_Constants/ Aux_const_ 1 104
```