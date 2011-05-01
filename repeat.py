#!/usr/bin/python

import os
import time
import sys

for x in range(5, 100):
   print('next')
   os.system('./make_input ./tests/dna.50MB ./tests/test.test 400 '+str(x)+' 0.2 -n 100000')
   if( os.system('./harness ./tests/test.test -periodic  > out.txt')  > 0):
      print('FAILED')
      break
   time.sleep(1)
