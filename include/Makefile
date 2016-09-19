# The point of this makefile is just to create the Main, and in the process
# update all the *.o files.

.KEEP_STATE:

OBJS_TIME = TimeMem.o
LIBS = libJtime.a libJgtools.a libJmarkov.a

# To use a library (e.g., libJtime) in a Makefile, do the following:
# INCLUDES = -I$(HOME)/include
# INC_LIBS = -L$(HOME)/include -lJtime
# then add $(INCLUDES) to the .cc.o dependencies and $(INC_LIBS) to the
# executable dependencies.

# compiler commands
CC      = g++
##CC = clang
RM      = /bin/rm -rf
AR      = ar crv
BACKUPS = *~ \\\#*\\\#

# compiler flags
CFLAGS  = -Wall -ansi -pedantic -g -O3 -std=c++11
INCLUDES = -I$(LACHESIS_SAMTOOLS_DIR) -I$(LACHESIS_BOOST_DIR)

# linking flags.  The LACHESIS_BOOST_DIR environment variable must be set
BOOST_LIBS=-L$(LACHESIS_BOOST_DIR)/stage/lib -lboost_filesystem
LFLAGS = $(BOOST_LIBS)

# dependencies

.cc.o:  .cc
	$(CC) -c $< $(CFLAGS) $(INCLUDES)

all:    $(LIBS)

libJtime.a: $(OBJS_TIME)
	$(AR) libJtime.a $(OBJS_TIME)

libJgtools.a:
	$(MAKE) -C gtools

libJmarkov.a:
	$(MAKE) -C markov

clean: 
	$(RM) $(LIBS) $(OBJS_TIME) gtools/*.o markov/*.o markov/TestMarkovModel core .make.state

clobber: clean
	$(RM) $(BACKUPS)

                                                                               
