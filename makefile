CC=g++
CFLAGS=-c -I.
LFLAGS=-g -lm -lgsl -lgslcblas 

%.out:%.o hdw.o conf.h
	$(CC) $(^:conf.h=) $(LFLAGS) -o $@

%.o:%.cpp hdw.conf
	$(CC) $(<:hdw.conf=) $(CFLAGS) -o $@

conf.h:hdw.conf
	./.conf2h $^ $@

clean:
	rm -rf *.o *.out
	rm -rf *~

PROGRAMS=hdw-test.out

all:$(PROGRAMS)

#######################################################################
#DEVELOPERS RULES 
#######################################################################

edit:
	emacs -nw *.cpp *.h *.conf MANUAL README

push:
	git commit -am "Commit"
	git push origin master

pull:
	git pull
