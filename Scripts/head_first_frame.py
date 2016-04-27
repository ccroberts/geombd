#!/usr/bin/env python

import sys, os

def tail(f, n):
  stdin,stdout = os.popen2('head -n %d %s' % (n, f))
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

for line in frame:
  sys.stdout.write(line)

