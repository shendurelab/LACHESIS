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

EXE = Lachesis
OBJS = Reporter.o ChromLinkMatrix.o GenomeLinkMatrix.o TrueMapping.o LinkSizeDistribution.o \
 ContigOrdering.o ClusterVec.o RunParams.o TextFileParsers.o Lachesis.o
CCFILES = Reporter.cc ChromLinkMatrix.cc GenomeLinkMatrix.cc TrueMapping.cc LinkSizeDistribution.cc \
 ContigOrdering.cc ClusterVec.cc RunParams.cc TextFileParsers.cc Lachesis.cc
BACKUPS = *~ \\\#*\\\#

Lachesis_CPPFLAGS = -I. -Iinclude $(SAMTOOLS_CPPFLAGS) $(BOOST_CPPFLAGS)
Lachesis_CFLAGS = -Wall -g -O3 -std=c++11 -ansi -pedantic -fPIC
Lachesis_LDFLAGS = -static -lz -lpthread -Linclude -lJtime -lJgtools -lJmarkov $(SAMTOOLS_LDFLAGS) $(LDFLAGS_BOOST) $(SAMTOOLS_LIBS)
Lachesis_SOURCES = $(CCFILES)
## I do not think I am supposed to need to add -lbam here, but that the macros from ax_lib_samtools.m4
## should provide it for me.
Lachesis_LDADD = $(AM_LDFLAGS) $(SAMTOOLS_LIBS) $(LIBS_BOOST) -lbam

bin_PROGRAMS = Lachesis
dist_bin_SCRIPTS = bin/CountMappables.pl bin/CountMotifsInFasta.pl \
 bin/CreateScaffoldedFasta.pl bin/PreprocessSAMs.pl bin/PreprocessSAMs.sh \
 bin/QuickDotplot bin/QuickDotplot.POA.R bin/QuickDotplot.R bin/QuickDotplot.SKY.R \
 bin/blast.qsub.sh bin/blast.sh bin/commentify_INIs.pl bin/decommentify_INIs.pl \
 bin/heatmap.MWAH.R bin/heatmap.R bin/make_bed_around_RE_site.pl bin/testme.sh \
 bin/INIs/test_case.ini

## This target may be used in order to invoke clang's tidy tool.
## It prints out a very nice report showing how one is or is not following the clang C++ style guide
## Use it with something like 'make ContigOrdering.tidy' and read the new file created by clang.
## Since I don't really know C++ at all, it provided a nice way for me to learn bits and pieces of the
## syntax which are not familar from C.
%.tidy: %.cc
	clang-tidy $< -checks=* -header-filter='.*' -- -Iinclude 2>$<.tidy 1>&2

LTest: $(OBJS) LTest.o
	$(CC) $(CFLAGS) $(OBJS) LTest.o -o LTest $(LFLAGS)
