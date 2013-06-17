#/usr/bin/env python

import sys
import random

f = open(sys.argv[1])
of = open(sys.argv[2],'w')
for l in f:
    nl = l.strip()+' '+str(random.randint(int(sys.argv[3]), int(sys.argv[4])))+'\n'
    of.write(nl)

f.close()
of.close()
