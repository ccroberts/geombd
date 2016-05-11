#!/usr/bin/python
import os, sys, struct

fileName = sys.argv[1]
fileContent = None

with open(fileName, mode='rb') as file:
  fileContent = file.read()

print struct.unpack("dddiiid", fileContent[:48])
