#!/bin/sh

# run a job
./exampleB4b -b paramBatch03_single.mac  -jobName testjob -runNumber 998 -runSeq 003 \
   -numberOfEvents 10  -eventsInNtuple 10 \
   -gun_particle e+ -gun_energy_min 100.0 -gun_energy_max 100.0 \
   -sipmType 1

