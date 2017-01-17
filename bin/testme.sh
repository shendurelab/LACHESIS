#!/usr/bin/env bash
export LACHESIS_SAMTOOLS_DIR=/cbcb/sw/RedHat-7-x86_64/common/local/samtools/0.1.19
export LACHESIS_BOOST_DIR=/cbcb/sw/RedHat-7-x86_64/common/local/boost/1.53
make clean
make -j5


./Lachesis INIs/test_case.ini 2>lach.err 1>lach.out
