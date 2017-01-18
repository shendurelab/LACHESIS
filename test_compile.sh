#!/usr/bin/env bash

echo "This silly script is intended to make it impossible for me to forget to re-invoke autoconf"

autoreconf -vif && ./configure --prefix=/tmp/test_lachesis_compile && make && make dist && make install && make uninstall && make distclean
