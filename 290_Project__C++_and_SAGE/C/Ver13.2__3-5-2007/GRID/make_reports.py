#!/usr/bin/python

import os
import time


def make_report(form_number):
    """ Read in an output file, and report the essential statsitics for the form computation. """

    # SANITY CHECKS
    assert form_number >= 1
    assert form_number <= 6560


    # Pre-declare all variables "defined" in conditional statements wiht nonsense values
    level = det = chi_top = 0
    aniso_vector_rawtext = ""
    basic_cusp_const = adjusted_cusp_const = -1

    # =========================================================================


    # Make the filename for this form
    form_filename = "OUTPUT_LOGS/Form" + str(form_number) + ".txt"


    # Read in the file as a string, and break it into lines
    form_file = open(form_filename, "r")
    output_lines = form_file.readlines()
    form_file.close()


    # Clean the lines of excess whitespace
    for i in range(len(output_lines)):
        output_lines[i] = output_lines[i].strip()
    

    # --------------------------------------------------------

    # Display the form number
    print "Summary for Form #" + str(form_number) + ":"
    print "-----------------------"
    print
    
    # --------------------------------------------------------

    print "Basic Invariants:"
    print "-----------------"

    # Look for the starting and ending time lines
    for line in output_lines:
        if line.startswith("Program started at:"):
            print line
        if line.startswith("Program finished at:"):
            print line

    # Print the Level (and Determinant and character)
    for line in output_lines:
        if line.startswith("The level is"):
            print line
            line = line.strip()
            level = int(line.lstrip("The level is "))
        if line.startswith("The determinant is"):
            print line
            line = line.strip()
            det = int(line.lstrip("The determinant is"))
        if line.startswith("The character is given by"):
            print line
            line = line.strip()
            chi_top = int(line.lstrip("The character is given by "))

    # Print the vector of anisotropic primes
    for line in output_lines:
        if line.startswith("The anisotropic primes are"):
            print line
            line = line.strip()
            aniso_vector_rawtext = line.lstrip("The anisotropic primes are ")
    print

    # --------------------------------------------------------

    print "Bounds:"
    print "-------"
    
    # Print the two cuspidal bounds
    for line in output_lines:
        if line.startswith("Passed in the cusp constant"):
            print line
            line = line.strip()
            basic_cusp_const = float(line.lstrip("Passed in the cusp constant "));
        if line.startswith("Using the larger cusp constant"):
            print line
            line = line.strip()
            adjusted_cusp_const = float(line.lstrip("Using the larger cusp constant ").rstrip(" to allow for roundoff errors! =)"))

    # Print the two Eisenstein bounds
    for line in output_lines:
        if line.startswith("The theoretical lower bound is:"):
            print line
            line = line.strip()
            theoretical_eis_const = 
        if line.startswith("The numerical lower bound is:"):
            print line
            line = line.strip()

    # Print the F4 bound
    for line in output_lines:
        if line.startswith("The F4 upper bound is:"):
            print line
            line = line.strip()

    print

    # --------------------------------------------------------

    print "Local Cover Info:"
    print "-----------------"

    # Print the local cover ternary
    for i in range(len(output_lines)):
        if output_lines[i].startswith("_form_list:"):
            for j in range(i+1, i+4):
                print output_lines[j]

    # Print the local cover unary coefficient
    for i in range(len(output_lines)):
        if output_lines[i].startswith("_d_list:"):
            print output_lines[i+1]

    print
    
    # --------------------------------------------------------

    print "Eligible Prime Info:"
    print "--------------------"

    # Print the number of eligible primes
    for line in output_lines:
        if line.startswith("There are") and line.endswith("eligible prime numbers."):
            print line
            line = line.strip()
    
    # Print the maximum number of eligible prime factors
    for line in output_lines:
        if line.startswith("We can have square-free numbers with at most") and line.endswith("prime factors."):
            print line
            line = line.strip()

    # Print the upper bound for eligible numbers
    for line in output_lines:
        if line.startswith("The largest eligible number must be less than"):
            print line
            line = line.strip()

    print

    # --------------------------------------------------------

    print "Number Checking Info:"
    print "---------------------"

    # Ternary precision used
    for line in output_lines:
        if line.startswith("using ternary precision"):
            print line
            line = line.strip()

    print

    # --------------------------------------------------------

    print "Square-free Exceptions:"
    print "-----------------------"

    # Print the square-free exceptions
    for line in output_lines:
        if line.startswith("The (square-free) exceptions of Q are:"):
            print line
            line = line.strip()
    
    print
    
    #============================================================


#    print
#    print "Python-extracted values:"
#    print "========================"
#
#    # Print the extracted variables
#    print "level =", level
#    print "det =", det
#    print "chi_top =", chi_top
#    print "aniso_vector_rawtext =", aniso_vector_rawtext
#
#    print "basic_cusp_const =", basic_cusp_const
#    print "adjusted_cusp_const =", adjusted_cusp_const
#
#    print "Exact (theoritcal) Eisenstein Lower bound =",
