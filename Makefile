.KEEP_STATE:

# source code
# EXES: individual binary executables
# OBJS: non-executable object files
EXES = Lachesis
OBJS = Reporter.o ChromLinkMatrix.o GenomeLinkMatrix.o TrueMapping.o LinkSizeDistribution.o ContigOrdering.o ClusterVec.o RunParams.o TextFileParsers.o
CCFILES = Reporter.cc ChromLinkMatrix.cc GenomeLinkMatrix.cc TrueMapping.cc LinkSizeDistribution.cc ContigOrdering.cc ClusterVec.cc RunParams.cc TextFileParsers.cc

# compiler commands
CC = g++
RM = /bin/rm -rf
BACKUPS = *~ \\\#*\\\#

# compiler flags
CFLAGS = -Wall -g -O3 -std=c++11
CFLAGS += -ansi
CFLAGS += -pedantic

# Includes
INCLUDE_DIR=include
INCLUDES=-I$(INCLUDE_DIR) -I$(LACHESIS_SAMTOOLS_DIR) -I$(LACHESIS_BOOST_DIR)

# Linking flags.  The environment variables BOOST_LIBS and SAMTOOLS_LIBS must be set or this won't work.
INC_LIBS=-L$(INCLUDE_DIR) -lJtime -lJgtools -lJmarkov
BOOST_LIBS=-L$(LACHESIS_BOOST_DIR)/stage/lib -lboost_system -lboost_filesystem -lboost_regex
SAMTOOLS_LIBS=-L$(LACHESIS_SAMTOOLS_DIR) -lbam
LFLAGS = $(INC_LIBS) $(BOOST_LIBS) $(SAMTOOLS_LIBS) -lz -lpthread

# dependencies
.cc.o:  .cc
	$(CC) -c $< $(CFLAGS) $(INCLUDES)

%.tidy: %.cc
	clang-tidy $< -checks=* -header-filter='.*' -- -Iinclude 2>$<.tidy 1>&2

all:   libs $(EXES)

libs:
	$(MAKE) -C include

Lachesis: $(OBJS) Lachesis.o
	$(CC) $(CFLAGS) $(OBJS) Lachesis.o -o Lachesis $(LFLAGS)

LTest: $(OBJS) LTest.o
	$(CC) $(CFLAGS) $(OBJS) LTest.o -o LTest $(LFLAGS)

clean:
	$(RM) $(EXES) $(OBJS) core .make.state gmon.out

clobber: clean
	$(RM) $(BACKUPS)
	(MAKE) clobber -C include # recurse to the include directory

# Environment variable check.
check-env:
ifndef LACHESIS_SAMTOOLS_DIR
    $(error Environment variable $$LACHESIS_SAMTOOLS_DIR is undefined - please set to a directory containing the samtools package)
endif
ifndef LACHESIS_BOOST_DIR
    $(error Environment variable $$LACHESIS_BOOST_DIR is undefined - please set to a directory containing the root of the Boost package)
endif
