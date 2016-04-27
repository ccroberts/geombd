#!/usr/bin/env python

import sys, os

def tail(f, n):
  stdin,stdout = os.popen2('tail -n %d %s' % (n, f))
  stdin.close()
  lines = stdout.readlines(); stdout.close()
  return lines

NlinesPerFrame = 0
Nligands = 0

for line in open(sys.argv[1], 'r'):
  NlinesPerFrame += 1
  if line.startswith('REMARK'):
    Nligands += 1
  if line.startswith('END'):
    break

frame = tail(sys.argv[1], NlinesPerFrame)

binder = False
for line in frame:
  if line.startswith('REMARK'):
    s = set(line.split())
    if 'bound=true' in s:
      binder = True
    else:
      if binder == True: print 'END'
      binder = False
  if binder == True:
    sys.stdout.write(line)

