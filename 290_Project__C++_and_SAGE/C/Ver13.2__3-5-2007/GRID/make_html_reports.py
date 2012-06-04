#!/usr/bin/python

import os
import time


def make_html_report(form_number):
    """ Read in an output file, and report the essential statsitics for the form computation. """

    # SANITY CHECKS
    assert form_number >= 1
    assert form_number <= 6560


    # Pre-declare all variables "defined" in conditional statements wiht nonsense values
    level = det = chi_top = 0
    aniso_vector_rawtext = ""
    basic_cusp_const = adjusted_cusp_const = -1

    # =========================================================================


    # Make the filenames for this form
    form_filename = "Basic_Output_Logs/Form" + str(form_number) + ".txt"
#    form_filename = "Auxiliary_Output_Logs/Form" + str(form_number) + ".txt"
#    form_filename = "OUTPUT_LOGS__1-1999/Form" + str(form_number) + ".txt"
    html_filename = "Basic_HTML_Reports/Form" + str(form_number) + ".html"
#    html_filename = "Auxiliary_HTML_Reports/Form" + str(form_number) + ".html"
#    html_filename = "HTML_REPORTS/Form" + str(form_number) + ".html"


    # Read in the file as a string, and break it into lines
    form_file = open(form_filename, "r")
    output_lines = form_file.readlines()
    form_file.close()


    # Open the output file for writing, and add some basics
    html_file = open(html_filename, "w")
    html_file.write("""<html><HEAD>

    <style type="text/css">
    <!--
    H1 {color: black;  text-align: center;  font-size: 24pt;  font-family: Arial}

    .HeadingText {font-size: 18pt;  font-family: Arial; font-weight: bold}
    .NormalText {text-indent: 2em; font-size: 12pt;  font-family: Arial;  font-style: italic}

    .SubMidItemText {font-size: 12pt;  font-family: Arial;  font-style: italic;
                     text-indent: 15px}
    .SubItemCorrectionText {font-size: 12pt;  font-family: Arial;  font-style: italic}
    .SubRefText {font-size: 12pt;  font-family: Arial;  font-style: italic;
                 text-indent: 25px}

    .MainFileText {font-size: 10pt;  font-family: Arial; text-indent: 25px}
    .SubFileText {font-size: 8pt;  font-family: Arial; text-indent: 25px}

    .RandomItemText {font-size: 16pt;  font-family: "Times New Roman"}
    -->
    </style>

    </HEAD> <body>""")
    html_file.write("<h1>Summary for Form #" + str(form_number) + ":</h1>")


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
    html_file.write('<div class="HeadingText">Basic Invariants</div>')

    # Look for the starting and ending time lines
    for line in output_lines:
        if line.startswith("Program started at:"):
            print line
            html_file.write('<div class="NormalText">' + line + '</div>')
        if line.startswith("Program finished at:"):
            print line
            html_file.write('<div class="NormalText">' + line + '</div>')
            
    # Print the Level (and Determinant and character)
    for line in output_lines:
        if line.startswith("The level is"):
            print line
            line = line.strip()
            level = int(line.lstrip("The level is "))
            html_file.write('<div class="NormalText">' + line + '</div>')
        if line.startswith("The determinant is"):
            print line
            line = line.strip()
            det = int(line.lstrip("The determinant is"))
            html_file.write('<div class="NormalText">' + line + '</div>')
        if line.startswith("The character is given by"):
            print line
            line = line.strip()
            chi_top = int(line.lstrip("The character is given by "))
            html_file.write('<div class="NormalText">' + line + '</div>')
            
    # Print the vector of anisotropic primes
    for line in output_lines:
        if line.startswith("The anisotropic primes are"):
            print line
            line = line.strip()
            aniso_vector_rawtext = line.lstrip("The anisotropic primes are ")
            html_file.write('<div class="NormalText">' + line + '</div>')
            
    print

    # --------------------------------------------------------

    print "Bounds:"
    print "-------"
    html_file.write('<div class="HeadingText">Relevant Bounds</div>')
    
    # Print the two cuspidal bounds
    for line in output_lines:
        if line.startswith("Passed in the cusp constant"):
            print line
            line = line.strip()
            basic_cusp_const = float(line.lstrip("Passed in the cusp constant "));
            html_file.write('<div class="NormalText">' + line + '</div>')
        if line.startswith("Using the larger cusp constant") and (basic_cusp_const != 0):
            print line
            line = line.strip()
            adjusted_cusp_const = float(line.lstrip("Using the larger cusp constant ").rstrip(" to allow for roundoff errors! =)"))
            html_file.write('<div class="NormalText">' + line + '</div>')
            
    # Print the two Eisenstein bounds
    for line in output_lines:
        if line.startswith("The theoretical lower bound is:"):
            print line
            line = line.strip()            
            html_file.write('<div class="NormalText">' + line + '</div>')
        if line.startswith("The numerical lower bound is:"):
            print line
            line = line.strip()
            html_file.write('<div class="NormalText">' + line + '</div>')
            
    # Print the F4 bound
    for line in output_lines:
        if line.startswith("The F4 upper bound is:"):
            print line
            line = line.strip()
            html_file.write('<div class="NormalText">' + line + '</div>')
            
    print

    # --------------------------------------------------------

    print "Local Cover Info:"
    print "-----------------"
    html_file.write('<div class="HeadingText">Local Cover Info</div>')

    # Print the local cover ternary
    for i in range(len(output_lines)):
        if output_lines[i].startswith("_form_list:"):
            for j in range(i+1, i+4):
                print output_lines[j]
                html_file.write('<div class="NormalText">' + output_lines[j] + '</div>')

    # Print the local cover unary coefficient
    for i in range(len(output_lines)):
        if output_lines[i].startswith("_d_list:"):
            print output_lines[i+1]
            html_file.write('<div class="NormalText">' + output_lines[i+1] + '</div>')

    print
    
    # --------------------------------------------------------

    print "Eligible Prime Info:"
    print "--------------------"
    html_file.write('<div class="HeadingText">Eligible Prime/Number Info</div>')

    # Print the number of eligible primes
    for line in output_lines:
        if line.startswith("There are") and line.endswith("eligible prime numbers."):
            print line
            line = line.strip()
            html_file.write('<div class="NormalText">' + line + '</div>')
                
    # Print the maximum number of eligible prime factors
    for line in output_lines:
        if line.startswith("We can have square-free numbers with at most") and line.endswith("prime factors."):
            print line
            line = line.strip()
            html_file.write('<div class="NormalText">' + line + '</div>')
            
    # Print the upper bound for eligible numbers
    for line in output_lines:
        if line.startswith("The largest eligible number must be less than"):
            print line
            line = line.strip()
            html_file.write('<div class="NormalText">' + line + '</div>')
            
    print

    # --------------------------------------------------------

    print "Number Checking Info:"
    print "---------------------"

    # Ternary precision used
    for line in output_lines:
        if line.startswith("using ternary precision"):
            print line
            line = line.strip()
            html_file.write('<div class="NormalText">' + line + '</div>')
            
    print

    # --------------------------------------------------------

    print "Square-free Exceptions:"
    print "-----------------------"
    html_file.write('<div class="HeadingText">Exceptions</div>')
    
    # Print the square-free exceptions
    for line in output_lines:
        if line.startswith("The (square-free) exceptions of Q are:"):
            print line
            line = line.strip()
            html_file.write('<div class="NormalText">' + line + '</div>')
                
    print

    # --------------------------------------------------------

    print "Timing Information:"
    print "-------------------"
    html_file.write('<div class="HeadingText">Timing Information</div>')

    # Look for the starting and ending time lines
    for line in output_lines:
        if line.startswith("Program started at:"):
            print line
            html_file.write('<div class="NormalText">' + line + '</div>')
        if line.startswith("Program finished at:"):
            print line
            html_file.write('<div class="NormalText">' + line + '</div>')
            


    
    #============================================================

    # Close the HTML file
    html_file.write('<a href="' + '../' + form_filename + '">raw output file</a></body></html>')
    html_file.close()


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


# Code to make HTML reports for a range of forms
#for i in range(1,2000):
#    make_html_report(i)       
