#!/usr/bin/env python

import sys, subprocess



def at_synonym(at):
  if at == 'H5\'': return 'H5\'1'
  elif at == 'H5\'\'': return 'H5\'2'
  elif at == 'H2\'': return 'H2\'1'
  elif at == 'H2\'\'': return 'H2\'2'
  elif at == 'HO5\'': return 'HO5\'1'
  elif at == 'HO5\'\'': return 'HO5\'2'
  elif at == 'HO2\'': return 'HO2\'1'
  elif at == 'HO2\'\'': return 'HO2\'2'
  return at


def at_to_element(at):
  return at[0]


def process_receptor(outf, af, ff, resname, resid, resdata):
  if resname in ff.keys():
    for line in resdata:
      at = line[12:16].strip()
      at = at_synonym(at)
      el = at_to_element(at)
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
      try:
        param = ff[resname][at]
        #print '%s  1.00  0.00   %8.4f %s' % (line[:54], param[0], param[1])
        pref = '%s%s%s' % (line[:12], str(' %s' % param[1].ljust(3, ' ')), line[16:30])
        outf.write('%s%10.4f%10.4f%10.4f %7.4f %6.4f\n' % (pref, x, y, z, param[0], af[param[1]]['r'])) #PQR
      except KeyError:
        print 'REMARK No charge assign for atom type', at
        pref = '%s%s%s' % (line[:12], str(' %s' % el.ljust(3, ' ')), line[16:30])
        #print '%s  1.00  0.00   %8.4f %s' % (line[:54], 0, at[0])
        outf.write('%s%10.4f%10.4f%10.4f %7.4f %6.4f\n' % (pref, x, y, z, 0, af[at[0]]['r'])) #PQR
  else:
    outf.write('REMARK No forcefield assignment for the following residue. Using OpenBABEL\n')
    tmpfd = open('/tmp/het.pdb', 'w')
    for line in resdata:
      tmpfd.write(line)
    tmpfd.close()
    subprocess.call(["babel", "-ipdb", "/tmp/het.pdb", "-opdb", "/tmp/heth.pdb"])
    subprocess.call(["babel", "-ipdb", "/tmp/heth.pdb", "-omol2", "/tmp/het.mol2"])
    Q = []
    tmpfd = open('/tmp/het.mol2', 'r')
    tmpfd.readline()
    tmpfd.readline()
    n = int(tmpfd.readline().split()[0])
    line = tmpfd.readline()
    while not line.startswith("@<TRIPOS>ATOM"):
      line = tmpfd.readline()
    for line in tmpfd:
      if line[0] == '@' or line.strip() == '': break
      sp = line.split()
      q = float(sp[-1])
      Q.append(q)
    i = 0
    for line in open('/tmp/heth.pdb', 'r'):
      if line.startswith('ATOM') or line.startswith('HETATM'):
        at = at_synonym(line[12:16].strip())
        el = at_to_element(at)
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        pref = '%s%s%s' % (line[:12], str(' %s' % el.ljust(3, ' ')), line[16:30])
        outf.write('%s%10.4f%10.4f%10.4f %7.4f %6.4f\n' % (pref, x, y, z, Q[i], af[el]['r'])) #PQR
        i += 1



def run(paramfn, pdbfn, outfn):
  ff = {}
  af = {}

  for line in open(paramfn, 'r'):
    sp = line.split()
    if len(sp) < 5: continue
    if sp[0] == 'frag_par':
      try:
        ff[sp[1]][sp[2]] = [float(sp[3]), sp[4]]
      except KeyError:
        ff[sp[1]] = {}
        ff[sp[1]][sp[2]] = [float(sp[3]), sp[4]]
    if sp[0] == 'atom_par':
      af[sp[1]] = {'r':float(sp[2])/2., 'eps':float(sp[3]), 'vol':float(sp[4]), 'sol':float(sp[5])}

  resn = None
  resi = -6969
  resd = []
  term = False

  outf = open(outfn, 'w')

  for line in open(pdbfn, 'r'):
    if line.startswith('ATOM') or line.startswith('HETATM'):
      rid = int(line[23:26])
      rnm = line[17:20].strip()
      if rid != resi:
        if resi == -6969:
          tnm = 'N%s' % rnm
          if tnm in ff.keys(): rnm = tnm
        else:
          process_receptor(outf, af, ff, resn, resi, resd)
          resd = []
          if term:
            outf.write('TER\n')
        if term:
          tnm = 'N%s' % rnm
          if tnm in ff.keys(): rnm = tnm
          term = False
        resi = rid
        resn = rnm
      resd.append(line)
      if line[12:16].strip() == 'OXT':
        resn = 'C%s' % resn
        term = True
    if line.startswith('END'):
      if len(resd) != 0:
        process_receptor(outf, af, ff, resn, resi, resd)
        resi = -6969
        resd = []
        resn = None
        outf.write('END\n')


if __name__ == '__main__':
  if len(sys.argv) < 7:
    print 'Usage:', sys.argv[0], '-d [Parm.gbdp]', '-i [Molecule.PDB]', '-o [Molecule.PQR]'
  if len(sys.argv) == 7:
    argd = { sys.argv[1]: sys.argv[2], sys.argv[3]: sys.argv[4], sys.argv[5]: sys.argv[6] }
    run(argd['-d'], argd['-i'], argd['-o'])
