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
RM      = /bin/rm -rf
AR      = ar crv
BACKUPS = *~ \\\#*\\\#

# directories with .h's in them
BOOST_DIR=/net/shendure/vol10/jnburton/extern/boost_1_47_0
INCLUDES = -I$(BOOST_DIR)

# compiler flags
CFLAGS  = -Wall -ansi -pedantic -g -O3

# linking flags
BOOST_LIBS=-L$(BOOST_DIR)/stage/lib -lboost_filesystem
LFLAGS = $(BOOST_LIBS)

# dependencies

.cc.o:  .cc
	$(CC) -c $< $(CFLAGS) $(INCLUDES)

all:    $(LIBS)

libJtime.a: $(OBJS_TIME)
	$(AR) libJtime.a $(OBJS_TIME)

libJgtools.a:
	cd gtools; make; cd ..

libJmarkov.a:
	cd markov; make; cd ..

clean: 
	$(RM) $(LIBS) $(OBJS_TIME) gtools/*.o markov/*.o markov/TestMarkovModel core .make.state

clobber: clean
	$(RM) $(BACKUPS)

                                                                               
