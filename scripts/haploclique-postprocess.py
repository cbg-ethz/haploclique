#! /usr/bin/env python
import sys
from cStringIO import StringIO
from math import pow, log10

prefix = ""
if len(sys.argv) == 2:
    prefix = "_"+sys.argv[1]

firstClique = False
buffer_lines = StringIO()

map = {}

def phredtodouble(s):
    letters = list(s) 
    letters = [pow(10,-1*float(ord(phred))/float(10)) for phred in letters] 
    return letters

def doubletophred(s):
    letters = StringIO()
    for x in s:
        c = int(-10*log10(float(x)/float(k.count))+33)
        if c > 126 or c < 33:
            c = 126
        letters.write(chr(c))
    return letters.getvalue()


def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)

def revcom(s):
    return complement(s[::-1])

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
    seq_1 = ""
    qual_1 = ""
    pos_1 = ""
    seq_2 = ""
    qual_2 = ""
    pos_2 = ""
    count = 1

key = ""
entry = Storage()
lineno = 0

SEP = "*|*"

for l in sys.stdin:
    if len(l.rstrip()) > 0:
        lineno += 1
        if l.startswith("@Clique1_"):
            if (len(entry.id) > 0):
                if firstClique:
                    key = entry.pos_1+SEP+entry.seq_1
                else:
                    key = entry.pos_1+SEP+entry.seq_1+SEP+entry.pos_2+SEP+entry.seq_2

                if key in map:
                    map[key].reads += entry.reads
                    map[key].reads = f5(map[key].reads)
                    map[key].count += 1
                    for x in xrange(len(map[key].qual_1)):
                        map[key].qual_1[x]+=entry.qual_1[x]
                    if len(map[key].qual_2) > 0:
                        for x in xrange(len(map[key].qual_2)):
                            map[key].qual_2[x]+=entry.qual_2[x]
                else:
                    map[key] = entry

            splits = l[:-1].split("-+-")
            entry = Storage()
            entry.reads = splits[1].split(",")
            entry.id = "Clique"+prefix+"_"+splits[0][9:]

            firstClique = True
        elif l.startswith("@Clique2_"):
            firstClique = False

        if not l.startswith("@Clique"):
            if firstClique:
                if lineno%4 == 2: entry.seq_1 = l[0:-1]
                if lineno%4 == 3: entry.pos_1 = l.split(" ")[0][1:]
                if lineno%4 == 0: entry.qual_1 = phredtodouble(l[0:-1])
            else:
                if lineno%4 == 2: entry.seq_2 = l[0:-1]
                if lineno%4 == 3: entry.pos_2 = l.split(" ")[0][1:]
                if lineno%4 == 0: entry.qual_2 = phredtodouble(l[0:-1])

if firstClique:
    key = entry.pos_1+SEP+entry.seq_1
else:
    key = entry.pos_1+SEP+entry.seq_1+SEP+entry.pos_2+SEP+entry.seq_2

if key in map:
    map[key].reads += entry.reads
    map[key].reads = f5(map[key].reads)
else:
    map[key] = entry

fr1 = open("data_cliques_paired_R1.fastq",'w')
fr2 = open("data_cliques_paired_R2.fastq",'w')
fs = open("data_cliques_single.fastq",'w')
fc = open("data_clique_to_reads.tsv", 'w')

for x in map:
    k = map[x]
    if k.pos_2 == "":
        fs.write("@")
        fs.write(k.id)
        fs.write("\n")
        fs.write(k.seq_1)
        fs.write("\n")
        fs.write("+")
        fs.write(k.pos_1)
        fs.write("\n")
        fs.write(doubletophred(k.qual_1))
        fs.write("\n")
    else:
        fr1.write("@")
        fr1.write(k.id)
        fr1.write("\n")
        fr1.write(k.seq_1)
        fr1.write("\n")
        fr1.write("+")
        fr1.write(k.pos_1)
        fr1.write("\n")
        fr1.write(doubletophred(k.qual_1))
        fr1.write("\n")
        fr2.write("@")
        fr2.write(k.id)
        fr2.write("\n")
        fr2.write(revcom(k.seq_2))
        fr2.write("\n")
        fr2.write("+")
        fr2.write(k.pos_2)
        fr2.write("\n")
        fr2.write(doubletophred(k.qual_2))
        fr2.write("\n")

    if len(k.reads) > 0:
        fc.write(k.id)
        fc.write("\t")
        fc.write(k.reads[0])
        for r in k.reads[1:]:
            fc.write(",")
            fc.write(r)
        fc.write("\n")

fr1.close()
fr2.close()
fs.close()
fc.close()
