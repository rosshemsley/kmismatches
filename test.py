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

# Check file does not already exist
if (os.path.exists(sys.argv[3])):
   print "File exists, continue? ",
   while (True):
      cont = sys.stdin.read(1)
      if (cont == 'n'): 
         sys.exit()
      if (cont == 'y'):
         break
   


# Flush the log before starting
#f= open('test.log', 'w')
#f.close()


for x in range(2,maxVal):
   this = []
   print 'Doing test ',x,' of ', maxVal
   for y in range(0,repeats):
      
      print 'Working on repeat ', y+1,' out of', repeats
      

      m = x*100
      k = x*5
      
      # Make a test case.
      
      maker = './make_input ./tests/english.50MB ./tests/test.test '+str(m)+' '+ str(k) + ' >>test.log'
      print maker
      os.system(maker)
        
        
      run   =   './km ./tests/test.test '+sys.argv[4]+' >>test.log'
      print run
      t0 = time.time();    
      # Run kmismatches
      os.system(run)
      
      this.append(time.time() - t0)

   average =  float(sum(this)) / len(this)
   
   stats = open(sys.argv[3], 'a')
   stats.write(str(x*100) + ' ' + str(average)+'\n')
   stats.close()
