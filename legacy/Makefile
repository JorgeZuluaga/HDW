CC = g++
#CFLAGS = -g -c -O4 -I.
CFLAGS = -g -c -I.
LFLAGS = -g -lm -lgslcblas -lgsl

clean:
	rm -rf *.o *.out *~ *\#

%.out:%.o
	$(CC) $^ $(LFLAGS) -o $@
