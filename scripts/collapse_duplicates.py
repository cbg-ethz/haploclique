#! /usr/bin/env python
import sys
from cStringIO import StringIO
from math import pow, log10

def phredtodouble(s):
    letters = list(s) 
    letters = [pow(10,-1*float(ord(phred))/float(10)) for phred in letters] 
    return letters

def f5(seq, idfun=None):
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

class Storage:
    id = ""
    reads = []
    line = []
    count = 1

key = ""
entry = Storage()
lineno = 0

SEP = "*|*"

map = {}

for l in sys.stdin:
    if len(l.rstrip()) > 0:
        entry = Storage()
        split = l.split(" ")
        if split[0].startswith("Clique"):
            first_split = split[0].split("-+-")
            entry.id = first_split[0]
            entry.reads = first_split[1].split(",")
        else:
            entry.id = split[0]
        entry.line = split[1:]

        key = split[6]+SEP+split[9]
        if len(split) > 17:
            key += SEP+split[14]+SEP+split[17]
        if key in map:
            map[key].reads += entry.reads
            map[key].reads = f5(map[key].reads)
            map[key].count += 1
        else:
            map[key] = entry


f = open("alignment.tmp",'w')

for x in map:
    k = map[x]
    f.write(k.id)
    if len(k.reads) > 0:
        f.write("-+-")
        f.write(",".join(k.reads))
    f.write(" ")
    f.write(" ".join(k.line))

f.close()
