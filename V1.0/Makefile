CFLAGS=-Wall -g -O3 -ffast-math  -march=native -fomit-frame-pointer -std=c99 
LFLAGS=-lgsl -lgslcblas -lm 
COMP=gcc
GGO=gengetopt

all : derrida         

hop1d_test : hop1d_test.c mt19937ar.o gaussian.o 
	$(COMP) hop1d_test.c mt19937ar.o gaussian.o $(CFLAGS) $(LFLAGS) -o hop1d_test

hop1d : hop1d.c mt19937ar.o gaussian.o 
	$(COMP) hop1d.c mt19937ar.o gaussian.o $(CFLAGS) $(LFLAGS) -o hop1d

hop1d_avg : hop1d_avg.c mt19937ar.o gaussian.o 
	$(COMP) hop1d_avg.c mt19937ar.o gaussian.o $(CFLAGS) $(LFLAGS) -o hop1d_avg

hop1d_avg_RB : hop1d_avg_RB.c mt19937ar.o gaussian.o 
	$(COMP) hop1d_avg_RB.c mt19937ar.o gaussian.o $(CFLAGS) $(LFLAGS) -o hop1d_avg_RB

hop1d_curves : hop1d_curves.c mt19937ar.o gaussian.o 
	$(COMP) hop1d_curves.c mt19937ar.o gaussian.o $(CFLAGS) $(LFLAGS) -o hop1d_curves

mt19937ar.o : mt19937ar.c
	$(COMP) -c mt19937ar.c -o mt19937ar.o $(CFLAGS)

gaussian.o : gaussian.c
	$(COMP) -c gaussian.c -ogaussian.o  $(CFLAGS)



clean :
	rm *.o hop1d
