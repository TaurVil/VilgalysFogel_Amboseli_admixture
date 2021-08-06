#!/bin/bash

# Concatenate all chromosome files per individual into a single file per individual
# h is the individual number
for h in `seq 1 508`; do cat $h.35kb.d2.* > $h.35kb.d2.masked.SWref.txt; done
