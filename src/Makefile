# Makefile for metsite
CC=gcc
CPP=g++
JAVAC=javac
LD=gcc
CFLAGS=-O3
LIBS=-lm
INC=/usr/include/

all: qmodcheck qmodcheck_mainens runBioserf tmjury_threadsafe

clean:
	rm -f qmodcheck
	rm -f qmodcheck_mainens
	rm -f runBioserf.class
	rm -f *subHitComparator.class
	rm -f ../bin/qmodcheck
	rm -f ../bin/qmodcheck_mainens
	rm -f ../bin/runBioserf*
	rm -f ../bin/tmjury3d_mq_modeller_threadsafe


install:
	cp qmodcheck ../bin
	cp qmodcheck_mainens ../bin
	cp runBioserf*.class ../bin
	cp tmjury3d_mq_modeller_threadsafe ../bin

qmodcheck:
		$(CC) $(CFLAGS) qmodcheck.c $(LIBS) -o qmodcheck

qmodcheck_mainens:
		$(CC) $(CFLAGS) qmodcheck_mainens.c $(LIBS) -o qmodcheck_mainens

runBioserf:
		$(JAVAC) -cp lib/biojava-1.7.1.jar:lib/bytecode.jar:./ runBioserf.java

tmjury_threadsafe:
		$(CC) $(CFLAGS) tmjury3d_mq_modeller_threadsafe.c $(LIBS) -o tmjury3d_mq_modeller_threadsafe
