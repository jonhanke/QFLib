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
##CXX=g++32
##CXX=g++
CXX=g++

LOCALDIR = ../../../Local
INCDIR = -I$(LOCALDIR)/include/pari -I$(LOCALDIR)/include -I/usr/local/include -I/.
LIBDIR = -L$(LOCALDIR)/lib -L/usr/local/lib
CFLAGS = $(INCDIR) $(LIBDIR) -O2 -lpari -lgmp -lgmpxx




## These are the C++ compiler flags  
CXXFLAGS= -I /usr/local/include -I /. \
	$(CFLAGS)

OBJ = distributed_class.o  misc_utilities.o  cholesky_distributed.o \
	Current_Local_Dir/Classes/Matrix_mpz.o \
	Current_Local_Dir/Classes/mpz_class_extras.o \
	boolean_theta_class.o \




DOXYGEN_DIR=Docs/Doxygen_Docs

#SAGE_DIR=~/public_html/SAGE/SAGE-qf-Devel__2005-05-30
#SAGE_DIR=~/public_html/SAGE/Current_Doxygen_Docs
WEBDOCS_DIR=~/public_html/QFLIB_Doxygen_Docs



all: main docs
## TO DO: Commented out for now -- Deal with this client code later (2013-12-09)
##all: main docs client client_new


Matrix_mpz_main: Matrix_mpz_main.o  Makefile  temp_clean  $(OBJ)
	$(CXX) $(CXXFLAGS) \
		-lm Matrix_mpz_main.o  -lgmp -lgmpxx -o Matrix_main $(OBJ)



#distributed_class.o: distributed_class.cc
#	$(CXX) $(CXXFLAGS) \
#	-lm distributed_class.cc  -o distributed_class.o:  



Boolean_theta_class_files = boolean_theta_class__defn.cc boolean_theta_class__methods.cc \
			boolean_theta_class__distributed_methods2.cc distributed_theta__client_more.cc \
			boolean_theta_class__distributed_methods_new.cc  multi-slice_client__more.cc




distributed_class.o:

misc_utilities.o:

cholesky_distributed.o: Current_Local_Dir/Classes/Matrix_mpz.o

Current_Local_Dir/Classes/Matrix_mpz.o:

Current_Local_Dir/Classes/mpz_class_extras.o:

boolean_theta_class.o: $(Boolean_theta_class_files)

#		$(CXX) $(CXXFLAGS) -c -o boolean_theta_class.o  boolean_theta_class.cc




# ----------------------------------------------------------------------------------------------------------------------------------



main__multi-slice_client.o: main__multi-slice_client.cc vec_utils.cc \
			boolean_theta_class.o 

main__theta-client2.o: XXX__Distributed_Stuff/main__theta-client2.cc vec_utils.cc \
			boolean_theta_class.o 



client: main__theta-client2.o Makefile
	$(CXX) $(CXXFLAGS) \
		-lm main__theta-client2.o boolean_theta_class.o -lgmp -lgmpxx -o ternary_theta_client -O2


client_new: main__multi-slice_client.o Makefile
	$(CXX) $(CXXFLAGS)  \
		-lm main__multi-slice_client.o boolean_theta_class.o  -lgmp -lgmpxx -o ternary_theta_client_new -O2




# ----------------------------------------------------------------------------------------------------------------------------------

main: clean temp_clean Makefile main.o
	$(CXX) $(CXXFLAGS)  \
		-lm main.o  -lpari -lgmp -lgmpxx -o main 



main_test: clean temp_clean Makefile main_test.o
	$(CXX) $(CXXFLAGS)  \
		-lm main_test.o  -lpari -lgmp -lgmpxx -o main_test 


# ----------------------------------------------------------------------------------------------------------------------------------




most:  main docs



temp_clean:
	rm -f ../TEMP_Distributed_Dir/*


docs: 
	#doxygen $(DOXYGEN_DIR)/Doxyfile
	rm -rf $(WEBDOCS_DIR)
	cp -rp $(DOXYGEN_DIR) $(WEBDOCS_DIR)


clean:
	rm -f *.o; rm -f *~; rm -f Current_Local_Dir/Classes/*.o

primes: main
	echo "Making the first 15 million primes...";
	echo "  Started at: ";
	date;
	./sieve -e 275604542 > primelist_upto_15millionth_prime.txt;
	echo "  Finished at: ";
	date	


