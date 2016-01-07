#!/usr/bin/env python

import sys

types = set([])

for line in open(sys.argv[1], 'r'):
  if line.startswith('ATOM') or line.startswith('HETATM'):
    sp = line.split()
    at = sp[-1]
    types.add(at)

print(types)
