OPTIONS = -O3 -Wall -march=native -ffast-math 
EXEC = tbd
INSTALLDIR := $(HOME)/bin/
COMPILER = gcc
# COMPILER = $(HOME)/bin/faster_clang

all: src/transportBuddy.c
	$(COMPILER) $(OPTIONS) -o $(EXEC) src/transportBuddy.c -lm
install : all
	cp $(EXEC) $(INSTALLDIR)
