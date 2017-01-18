# The point of this makefile is just to create the Main, and in the process
# update all the *.o files.

.KEEP_STATE:

OBJS = N50.o ChromInterval.o HumanGenome.o CopyNumberProfile.o FileParsers.o SAMStepper.o
LIB  = libJgtools.a

# compiler commands
CC      = g++
##CC = clang
RM      = /bin/rm -rf
AR      = ar crv
BACKUPS = *~ \\\#*\\\#

# directories with .h's in them - the LACHESIS_BOOST_DIR and LACHESIS_SAMTOOLS_DIR environment variables must be set
INCLUDES = 

# compiler flags
CFLAGS  = -Wall -ansi -pedantic -g -O3 -std=c++11

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

                                                                               
