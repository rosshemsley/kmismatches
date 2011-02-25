#!/usr/bin/python

import os
import sys
import time

if(len(sys.argv) < 3):
   print 'Too few arguments. Exiting.\n'
   sys.exit()


repeats = int(sys.argv[1])
maxVal  = int(sys.argv[2])

print 'Writing output to "', sys.argv[3], '"\n'

# Clear the contents of the output file.
f = open(sys.argv[3], 'w');
f.close()


for x in range(2,maxVal):
   this = []
   for y in range(0,repeats):
   
      t0 = time.time();

      # Run kmismatches
      os.system('./km ./tests/case2/test'+str(x)+'.test '+sys.argv[3]+' >>test.log')
      
      this.append(time.time() - t0)

   average =  float(sum(this)) / len(this)
   
   stats = open(sys.argv[3], 'a')
   stats.write(str(x*100) + ' ' + str(average)+'\n')
   stats.close()
