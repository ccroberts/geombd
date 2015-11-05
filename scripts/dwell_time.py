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

NnoContact = 0
Ncontact = 0
Nbind = 0
tAvgDwell = 0.

for line in frame:
  if line.startswith('REMARK'):
    params = [s.strip().split('=') for s in line.split()[1:]]
    if len(params) == 0: continue
    if params[0][0] == 'done' and len(params) == 1:
      NnoContact += 1
    if len(params) > 1:
      Ncontact += 1
      tAvgDwell += float(params[1][1])
    if params[0][0] == 'bound':
      Nbind += 1

print('divided by N_contacts:', (tAvgDwell / Ncontact))
print('divided by N_ligands:', (tAvgDwell / Nligands))
print('N_contacts:', Ncontact)
print('N_no_contacts:', NnoContact)
print('N_bind:', Nbind)
print('tAvgDwell:', tAvgDwell)
