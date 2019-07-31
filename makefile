OPTIONS = -O2 -Wall
EXEC = tbd
INSTALLDIR := $(HOME)/bin/
COMPILER = gcc

all: src/transportBuddy.c
	$(COMPILER) $(OPTIONS) -o $(EXEC) src/transportBuddy.c
install : all
	cp $(EXEC) $(INSTALLDIR)