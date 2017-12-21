### Notes for the installation on GMI HPC cluster MENDEL

the implemented fixes should allow the simple ./configure && make && make install

# used dependencies
# ... using an older version of samtools (0.1.19 or earlier)

ml Autotools/20150215-foss-2016a
ml SAMtools/0.1.19-foss-2016a
ml Boost/1.60.0-foss-2016a

export LACHESIS_BOOST_DIR=$EBROOTBOOST
export LACHESIS_SAMTOOLS_DIR=$EBROOTSAMTOOLS

./configure --with-samtools=$EBROOTSAMTOOLS --with-boost=$EBROOTBOOST

# === THE END ===

#####################################
# the steps below were the reasoning how to reach this point
#####################################

# ./configure script fails with this error:
# configure: error: cannot run /bin/sh ./config.sub

# re-running autotools


# install missing autoconf files
autoreconf -i

# now configure gets further, but fails with:
# checking samtools... failed
# configure: error: either specify a valid samtools installation with --with-samtools=DIR or disable samtools usage with --without-samtools

# funny sidenote: the LACHESIS_*_DIR variable are only found in README and ChangeLog, bot not in the configure script. need to use the --with-*=DIR flags

# ./configure --help shows a variable BOOST_ROOT
# is it expecting Boost 1.52.0 ?
# --with-boost=DIR        prefix of Boost 1.52.0 [guess]

./configure --with-samtools=$EBROOTSAMTOOLS --with-boost=$EBROOTBOOST

# there is also a switch for linking to static boost libraries... (but we did not use it)

# configure works,
# do 
find ./ -iname '*.a' -type f -delete 
find ./ -iname '*.o' -type f -delete 

# to remove all remaining object files, re-run make
# in SAMStepper.cc SAMStepper.h fix include of sam.h to bam/sam.h

gives some cx11 compile error...

# fix CFLAGS in src/include/markov/Makefile
# from CFLAGS = -Wall -ansi -pedantic -g -O3 -std=c++11
# to CFLAGS = -Wall -ansi -pedantic -g -O3

# it seems to compile then
