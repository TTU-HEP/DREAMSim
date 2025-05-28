#!/bin/bash

c++ partonShower.cc -o test `root-config --cflags --libs`   -I$PYTHIA8/include -L$PYTHIA8/lib -lpythia8

