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

# Includes
INCLUDE_DIR=include
INCLUDES=-I$(INCLUDE_DIR) -I$(LACHESIS_SAMTOOLS_DIR)

# Linking flags.  The environment variables BOOST_LIBS and SAMTOOLS_LIBS must be set or this won't work.
INC_LIBS=-L$(INCLUDE_DIR) -lJtime -lJgtools -lJmarkov
BOOST_LIBS=-L$(LACHESIS_BOOST_DIR)/stage/lib -lboost_system -lboost_filesystem -lboost_regex
SAMTOOLS_LIBS=-L$(LACHESIS_SAMTOOLS_DIR) -lbam
LFLAGS = $(INC_LIBS) $(BOOST_LIBS) $(SAMTOOLS_LIBS) -lz -lpthread


# dependencies

.cc.o:  .cc
	$(CC) -c $< $(CFLAGS) $(INCLUDES)

all:   check-env libs $(EXES)

libs: check-env
	$(MAKE) -C include

Lachesis: check-env $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o Lachesis $(LFLAGS)

clean:
	$(RM) $(OBJS) core .make.state gmon.out
	cd include; make clean; cd ..

clobber: clean
	$(RM) $(BACKUPS) $(EXES)


# Environment variable check.
check-env:
ifndef LACHESIS_SAMTOOLS_DIR
    $(error Environment variable $$LACHESIS_SAMTOOLS_DIR is undefined - please set to a directory containing the samtools package)
endif
ifndef LACHESIS_BOOST_DIR
    $(error Environment variable $$LACHESIS_BOOST_DIR is undefined - please set to a directory containing the root of the Boost package)
endif

                                                                               
