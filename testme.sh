#!/usr/bin/env bash
cd $(dirname $(which Lachesis))

./Lachesis INIs/test_case.ini 2>lach.err 1>lach.out
