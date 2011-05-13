#!/usr/bin/env python

import os
import sys
import time

if(len(sys.argv) < 8):
   print 'Expected [data, start, end, interval repeats, k, outfile, args]\n'
   sys.exit()

data     =       sys.argv[1]
start    =   int(sys.argv[2])
end      =   int(sys.argv[3])
interval =   int(sys.argv[4])
repeats  =   int(sys.argv[5])
k_ratio  = float(sys.argv[6])
output   =       sys.argv[7]

# If there is an argument to pass to the program.  
if (len(sys.argv)>8):
   arg   =       sys.argv[8]
else:
   arg = ''

log = output + '.log'

print 'Writing output to', output

## Check file does not already exist
#if (fileexistsoutput):
#   print "File exists, continue? ",
#   while (True):
#      cont = sys.stdin.read(1)
#      if (cont == 'n'): 
#         sys.exit()
#      if (cont == 'y'):
#         break



# Make and run the tests using the following.
run   = './harness .temp.test '+arg+' >>'+log

# Flush the output files before starting.
f= open(output, 'w')
f.close()

f= open(log, 'w')
f.write('# Original Data:'+ data+ ' Args: '+arg+'\n');
f.close()

# Run the tests.
for x in range(start,end,interval):
   
   # The array where we store the results for each run.
   this = []
   
   print 'Doing test ',(x-start+interval)/interval,' of ', (end-start)/interval
   
   for y in range(0,repeats):
      
      print 'Working on repeat ', y+1,' out of', repeats
      
      m = x
      k = x*k_ratio
      
      maker = './make_input '+data+' .temp.test '+str(m)+' '+ str(k) +' 0.01  >>'+log
       
      print 
      print ' $', maker

      if (os.WEXITSTATUS(os.system(maker)) != 0):
         print 'Failed to construct a test'
         exit(1)          

      # Give some time for system resources to be returned etc.
      time.sleep(2)
                        
      print ' $', run 
      
      t0 = time.time();
      
      if (os.WEXITSTATUS(os.system(run)) != 0):
         print 'A Test Failed'
         exit(1)    


      print               
      this.append(time.time() - t0)
      

     
   average =  float(sum(this)) / len(this)
   
   stats = open(output, 'a')
   stats.write(str(x) + ' ' + str(average) + ' ')
   
   # Write out all the values we found after the average.
   for v in this:
      stats.write( str(v) + ' ')
   stats.write('\n');
   
   
    
   stats.close()
   os.system('rm ./.temp.test')
   
   
