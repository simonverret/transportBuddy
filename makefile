OPTIONS = -O2 -Wall
EXEC = tbd
INSTALLDIR := $(HOME)/bin/
#COMPILER = gcc
COMPILER = /Users/Simon/codes/faster_clang

all: executable
executable: transportBuddy.c
	$(COMPILER) $(OPTIONS) -o $(EXEC) transportBuddy.c -llapack -lblas -lstdc++ -lm
install : all
	cp $(EXEC) $(INSTALLDIR)