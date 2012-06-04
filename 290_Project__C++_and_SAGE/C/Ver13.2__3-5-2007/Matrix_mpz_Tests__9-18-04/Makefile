######################################
#				     #
# This is the makefile for Version 8 #
# of my Quadratic Forms project.     #
#				     #
######################################
#			#
# Cleaned up on 9-9-04	#
#			#
#########################


## This is unnecessary because of the implicit rules
CXX=g++

## These are the C++ compiler flags  
CXXFLAGS= -I /usr/local/include -I /. -L /usr/local/lib  -lgmp -lgmpxx -g
#CXXFLAGS= -Wall

OBJ = Matrix_mpz.o LCM.o






all: Matrix_mpz_main


Matrix_mpz_main: Matrix_mpz_main.o  Makefile  $(OBJ)
	g++ $(CXXFLAGS) -I /usr/local/include -I /. -L /usr/local/lib  \
		-lm Matrix_mpz_main.o  -lgmp -lgmpxx -o Matrix_main $(OBJ)

Matrix_mpz.o: Matrix_mpz.cc

LCM.o: LCM.cc

clean:
	rm -f *.o

