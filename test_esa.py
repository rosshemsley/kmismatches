#!/usr/bin/env python

import os, time

this = []
for i in range (0, 4):
   t0 = time.time();    
   if (os.WEXITSTATUS(os.system('./test')) != 0):
      print 'A Test Failed'
      exit(1)
   this.append(time.time() - t0)   
   
average =  float(sum(this)) / len(this)   
print this
print average
