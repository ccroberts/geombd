import sys

for line in open(sys.argv[1]):
  sp = line.split()
  print '%8s %8s %16.6f %8s' % (sp[0], sp[1], float(sp[2]), sp[3])
