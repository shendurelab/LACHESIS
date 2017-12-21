.KEEP_STATE:

# source code
# EXES: individual binary executables
# OBJS: *all* non-executable object files (all EXEs require all OBJs; I don't bother to unravel the dependency web)
EXES = TestMarkovModel
OBJS = SymbolSet.o WDAG.o MarkovModel.o MarkovChain.o HMM.o
LIB  = libJmarkov.a

# binary (executable) path

# compiler commands

CC      = g++
##CC = clang
RM      = /bin/rm -rf
AR      = ar crv
BACKUPS = *~ \\\#*\\\#

# directories with .h's in them - the LACHESIS_BOOST_DIR environment variable must be set
INCLUDES =

# compiler flags

CFLAGS  = -Wall -ansi -pedantic -g -O3 
#CFLAGS += -pg

# linking flags
BOOST_LIBS=-lboost_system -lboost_filesystem -lboost_regex
#INC_LIBS=-L$(HOME)/include -lJtime
LFLAGS = -lz $(BOOST_LIBS) $(INC_LIBS)

# dependencies

.cc.o:  .cc
	$(CC) -c $< $(CFLAGS) $(INCLUDES)

all:    $(EXES) $(LIB)

TestMarkovModel:  TestMarkovModel.o $(OBJS)
	$(CC) $(CFLAGS) $< $(OBJS) -o TestMarkovModel $(LFLAGS)

$(LIB): $(OBJS)
	$(AR) $(LIB) $(OBJS); mv $(LIB) ..

clean: 
	$(RM) $(OBJS) core .make.state

clobber: clean
	$(RM) $(BACKUPS) $(EXES)

                                                                               
