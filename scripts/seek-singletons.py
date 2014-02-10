#!/usr/bin/env python
import sys

#get fname from parameter
idfile=sys.argv[1]

#load ids
ids = []
for x in open(idfile):
  ids += x.strip().split("\t")[1].split(",")
ids = set(ids)

for line in sys.stdin:
  if line.startswith("Clique"): continue
  id=line.split(" ")[0]
  if not id in ids:
    print(line.rstrip())
