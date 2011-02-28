#!/usr/bin/python

import os
import sys
import time

if(len(sys.argv) < 4):
   print 'Too few arguments. Exiting.\n'
   sys.exit()

repeats = int(sys.argv[1])
maxVal  = int(sys.argv[2])
data    = './tests/english.50MB'

print 'Writing output to', sys.argv[3]

# Check file does not already exist
if (os.path.exists(sys.argv[3])):
   print "File exists, continue? ",
   while (True):

      cont = sys.stdin.read(1)
      if (cont == 'n'): 
         sys.exit()
      if (cont == 'y'):
         break

# If there is an argument to pass to the program.  
if (len(sys.argv)>4):
   arg = sys.argv[4]
else:
   arg = ''

# Make and run the tests using the following.
run   = './km ./tests/test.test '+arg+' >>test.log'

# Flush the output files before starting.
f= open('test.log', 'w')
f.close()

f= open(sys.argv[3], 'w')
f.write('# Original Data:'+ data+ ' Args: '+arg+'\n');
f.close()

# Run the tests.
for x in range(2,maxVal+1):
   
   # The array where we store the results for each run.
   this = []
   
   print 'Doing test ',x,' of ', maxVal
   
   for y in range(0,repeats):
      
      print 'Working on repeat ', y+1,' out of', repeats
      
      m = x*100
      k = x*5
      
      maker = './make_input '+data+' ./tests/test.test '+str(m)+' '+ str(k) + ' >>test.log'
       
      print 
      print ' $', maker
      os.system(maker)
             
      print ' $', run
      print 
      
      t0 = time.time();    
      os.system(run)      
      this.append(time.time() - t0)

     
   average =  float(sum(this)) / len(this)
   
   stats = open(sys.argv[3], 'a')
   stats.write(str(x*100) + ' ' + str(average) + ' ')
   
   # Write out all the values we found after the average.
   for v in this:
      stats.write( str(v) + ' ')
   stats.write('\n');
   
   
    
   stats.close()
   
   
   
