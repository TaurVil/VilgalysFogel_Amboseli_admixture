#!/bin/bash

# For chromosomes 1-20 (defined by g)
# NUMBER defines individual
for g in `seq 1 20`; do sed -i 's/^/\t/' NUMBER.35kb.d2.$g.*txt; sed -i s/^/$g/g NUMBER.35kb.d2.$g.masked.SWref.txt; done
