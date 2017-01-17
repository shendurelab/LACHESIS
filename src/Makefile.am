ACLOCAL_AMFLAGS = -I m4 --install

## These flags are being extracted via macros found in m4/boost.m4
LDFLAGS_BOOST = \
    $(BOOST_SYSTEM_LDFLAGS) \
    $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
## Ditto, the variables created herein will be used by the variables Lachesis_* in order to
## properly include the libraries used from boost.
LIBS_BOOST = \
    $(BOOST_SYSTEM_LIBS) \
    $(BOOST_PROGRAM_OPTIONS_LIBS) \
    $(BOOST_FILESYSTEM_LIBS)

SAMTOOLS_CPPFLAGS = -Iinclude/samtools
SAMTOOLS_LDFLAGS = -Linclude/samtools -lbam
SAMTOOLS_LIBS = -lbam

EXE = Lachesis
OBJS = Reporter.o ChromLinkMatrix.o GenomeLinkMatrix.o TrueMapping.o LinkSizeDistribution.o \
 ContigOrdering.o ClusterVec.o RunParams.o TextFileParsers.o Lachesis.o
CCFILES = Reporter.cc ChromLinkMatrix.cc GenomeLinkMatrix.cc TrueMapping.cc LinkSizeDistribution.cc \
 ContigOrdering.cc ClusterVec.cc RunParams.cc TextFileParsers.cc Lachesis.cc
BACKUPS = *~ \\\#*\\\#

bin_PROGRAMS = Lachesis
Lachesis_CPPFLAGS = -I. -Iinclude $(SAMTOOLS_CPPFLAGS) $(BOOST_CPPFLAGS)
Lachesis_CFLAGS = -Wall -g -O3 -std=c++11 -ansi -pedantic -fPIC
Lachesis_LDFLAGS = -static -lz -lpthread -Linclude -lJtime -lJgtools -lJmarkov $(SAMTOOLS_LDFLAGS) $(LDFLAGS_BOOST)
Lachesis_SOURCES = $(CCFILES)
Lachesis_LDADD = $(AM_LDFLAGS) $(SAMTOOLS_LIBS) $(LIBS_BOOST)

## dependencies
##.cc.o:  .cc
##	$(CC) -c $< $(CFLAGS) $(INCLUDES)
##
## This target may be used in order to invoke clang's tidy tool.
## It prints out a very nice report showing how one is or is not following the clang C++ style guide
## Use it with something like 'make ContigOrdering.tidy' and read the new file created by clang.
## Since I don't really know C++ at all, it provided a nice way for me to learn bits and pieces of the
## syntax which are not familar from C.
##%.tidy: %.cc
##	clang-tidy $< -checks=* -header-filter='.*' -- -Iinclude 2>$<.tidy 1>&2
##
##all:   libs $(EXES)
##
##libs:
##	$(MAKE) -C include
##
##Lachesis: $(OBJS) Lachesis.o
##	$(CC) $(CFLAGS) $(OBJS) Lachesis.o -o Lachesis $(LFLAGS)
##
##LTest: $(OBJS) LTest.o
##	$(CC) $(CFLAGS) $(OBJS) LTest.o -o LTest $(LFLAGS)
##
##clean:
##	$(RM) $(EXES) $(OBJS) core .make.state gmon.out
##
##clobber: clean
##	$(RM) $(BACKUPS)
##	(MAKE) clobber -C include # recurse to the include directory
