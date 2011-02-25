#!/usr/bin/python

import os
import time



x = 2

for x in range(2,3):

   # Create a new test case.
   os.system('./make_input ./tests/english.50MB ./tests/case2/test'+str(x)+'.test ' + str(x*100) + ' ' + str(x*10))


   log = open('./test.log', 'a')
   log.write('About to do test '+str(x)+'\n')
   log.close();

   t0 = time.time();

   # Run kmismatches
   os.system('./km ./tests/case2/test'+str(x)+'.test 2>>test.log')

   km = time.time() - t0
   t0 = time.time();

   # Run naive
   os.system('./km ./tests/case2/test'+str(x)+'.test -naive 2>>test.log')
   naive = time.time() - t0

   print x, ': KM took ', km, ' Naive took ', naive

   stats = open('./stats', 'a')
   stats.write(str(x*100) + ' ' + str(km) + ' ' + str(naive) + '\n')
   stats.close()
