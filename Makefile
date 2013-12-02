.KEEP_STATE:

# source code
# EXES: individual binary executables
# OBJS: non-executable object files
EXES = Lachesis
OBJS = Lachesis.o Reporter.o ChromLinkMatrix.o GenomeLinkMatrix.o TrueMapping.o LinkSizeDistribution.o ContigOrdering.o ClusterVec.o RunParams.o TextFileParsers.o

# compiler commands

CC      = g++
RM      = /bin/rm -rf
BACKUPS = *~ \\\#*\\\#


# compiler flags
CFLAGS = -Wall -g -O3
CFLAGS += -ansi
CFLAGS += -pedantic
#CFLAGS += -pg

# Directories for includes & linking flags
INCLUDE_DIR=include
BOOST_DIR=/net/shendure/vol10/jnburton/extern/boost_1_47_0
SAMTOOLS_DIR=/net/shendure/vol10/jnburton/extern/samtools-0.1.16

# Includes
INCLUDES=-I$(INCLUDE_DIR) -I$(SAMTOOLS_DIR)

# Linking flags
INC_LIBS=-L$(INCLUDE_DIR) -lJtime -lJgtools -lJmarkov
BOOST_LIBS=-L$(BOOST_DIR)/stage/lib -lboost_system -lboost_filesystem -lboost_regex
SAMTOOLS_LIBS=-L$(SAMTOOLS_DIR) -lbam
LFLAGS = $(INC_LIBS) $(BOOST_LIBS) $(SAMTOOLS_LIBS) -lz


# dependencies

.cc.o:  .cc
	$(CC) -c $< $(CFLAGS) $(INCLUDES)

all:    libs $(EXES)

libs:
	cd $(INCLUDE_DIR); make; cd ..

Lachesis:  $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o Lachesis $(LFLAGS)

clean:
	$(RM) $(OBJS) core .make.state gmon.out
	cd include; make clean; cd ..

clobber: clean
	$(RM) $(BACKUPS) $(EXES)

                                                                               
