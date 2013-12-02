# The point of this makefile is just to create the Main, and in the process
# update all the *.o files.

.KEEP_STATE:

OBJS = N50.o ChromInterval.o HumanGenome.o CopyNumberProfile.o FileParsers.o SAMStepper.o
LIB  = libJgtools.a



# compiler commands
CC      = g++
RM      = /bin/rm -rf
AR      = ar crv
BACKUPS = *~ \\\#*\\\#

# directories with .h's in them
BOOST_DIR=/net/shendure/vol10/jnburton/extern/boost_1_47_0
SAMTOOLS_DIR=/net/shendure/vol10/jnburton/extern/samtools-0.1.16
INCLUDES = -I$(BOOST_DIR) -I$(SAMTOOLS_DIR)

# compiler flags
CFLAGS  = -Wall -ansi -pedantic -g -O3

# linking flags - these don't do anything because there's no linking in compilation
#BOOST_LIBS=-L$(BOOST_DIR)/stage/lib -lboost_filesystem -lboost_regex
#LFLAGS = $(BOOST_LIBS)

# dependencies

.cc.o:  .cc
	$(CC) -c $< $(CFLAGS) $(INCLUDES)

all:    $(LIB)

$(LIB): $(OBJS)
	$(AR) $(LIB) $(OBJS); mv $(LIB) ..

clean: 
	$(RM) $(LIBS) $(OBJS) core .make.state

clobber: clean
	$(RM) $(BACKUPS)

                                                                               
