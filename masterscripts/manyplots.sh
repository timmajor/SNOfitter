#!/bin/bash

for i in {0..0}
do
  python plotFit.py "$i" 0
  python plotFit.py "$i" 2 "all"
done
#python plotFit.py 1 2 "all"
#python plotFit.py 2 2 "all"
