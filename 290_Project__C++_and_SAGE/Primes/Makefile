CC=gcc

CXXFLAGS=-Wall

all: main

main.o: sieve2310_64bit.c

main: sieve2310_64bit.o
	gcc $(CXXFLAGS) -I /usr/local/include -I /. -L /usr/local/lib  \
		-lm sieve2310_64bit.o -o sieve -O2

primes15: main
	echo "Making the first 15 million primes...";
	echo "  Started at: ";
	date;
	./sieve -e 275604542 > primelist_upto_15millionth_prime.txt;
	echo "  Finished at: ";
	date	

primes40: main
	echo "Making the first 40 million primes...";
	echo "  Started at: ";
	date;
	./sieve -e 776531402 > primelist_upto_40millionth_prime.txt;
	echo "  Finished at: ";
	date	
