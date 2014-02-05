#! /usr/bin/env python
import sys
from cStringIO import StringIO

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

map = {}
id = ""
seq = ""
qual = ""
pos = ""
lineno = 0

singles = open('data_cliques_single.fastq', 'r')
for l in singles:
    lineno += 1
    if lineno%4 == 1:
    	if len(seq) > 0:
            key = pos+"-"+seq
            if key in map:
                map[key] += "*|*"+id
            else:
                map[key] = id
        id = l[0:-1]
    if lineno%4 == 2: seq = l[0:-1]
    if lineno%4 == 3:
    	pos = l.split(" ")[0][1:]
    if lineno%4 == 0: qual = l

if len(seq) > 0:
    key = pos+"-"+seq
    if key in map:
        map[key] += "*|*"+id
    else:
        map[key] = id

file = open('data_clique_to_reads.tsv', 'r')
tab = {}
for l in file.readlines():
    split = l.split("\t")
    tab[split[0]] = split[1][:-1]

fastq = StringIO()
ids = StringIO()
for key in map:
    split = key.split("-")
    idSplit = map[key].split("*|*")
    fastq.write(idSplit[0]+"\n")
    fastq.write(split[1]+"\n")
    fastq.write("+"+split[0]+"\n")
    for l in xrange(len(split[1])):
        fastq.write("~")
    fastq.write("\n")
    duplicate_ids = []
    for x in idSplit:
        for y in tab[x[1:]].split(","):
            duplicate_ids.append(y)
    duplicate_ids = f5(duplicate_ids)
    ids.write(idSplit[0][1:]+"\t"+duplicate_ids[0])
    tab.pop(idSplit[0][1:],None)
    for x in idSplit[1:]:
        tab.pop(x[1:],None)
    for x in duplicate_ids[1:]:
        ids.write(","+x)
    ids.write("\n")

for x in tab:
    ids.write(x+"\t"+tab[x]+"\n")

f = open("data_clique_to_reads.tsv", 'w')
f.write(ids.getvalue())
f.close()

f = open("data_cliques_single.fastq",'w')
f.write(fastq.getvalue())
f.close()
