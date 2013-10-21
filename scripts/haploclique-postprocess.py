#! /usr/bin/env python
import sys
from cStringIO import StringIO

sequences_paired1 = StringIO()
sequences_paired2 = StringIO()
sequences_single = StringIO()

firstClique = False
fastq = 0

out = StringIO()
name = ""
seq = ""
found = 0
buffer_lines = StringIO()

readnames = {}
countlist = []

clique_to_reads = StringIO()

for lines in sys.stdin:
    if len(lines) > 0:
        if lines.startswith("@Clique1_"):
            if firstClique:
                sequences_single.write(buffer_lines.getvalue())
                buffer_lines = StringIO()
            else:
                sequences_paired2.write(buffer_lines.getvalue())
                buffer_lines = StringIO()
            splits = lines[1:-1].split("-+-")
            clique_to_reads.write(splits[0].replace("Clique1","Clique"))
            clique_to_reads.write("\t");

            readSplit = splits[1].split(",")
            readCount = {}
            for r in readSplit:
                name_tmp = ""
                count_tmp = 1
                if "-" in r:
                    split_tmp = r.split("-")
                    name_tmp = split_tmp[1]
                    count_tmp = split_tmp[0]
                else:
                    name_tmp = r
                if name_tmp in readCount:
                    readCount[name_tmp] += count_tmp
                else:
                    readCount[name_tmp] = count_tmp
                if name_tmp in readnames:
                    readnames[name_tmp] += count_tmp
                else:
                    readnames[name_tmp] = count_tmp
            start = 1
            for r in readCount:
                if r != "1" and r != "":
                    if (start == 0):
                        clique_to_reads.write(",")
                    else:
                        start = 0
                    clique_to_reads.write(str(readCount[r]))
                    clique_to_reads.write("-")
                    clique_to_reads.write(r)
            clique_to_reads.write("\n")

            countlist.append(readCount)

            name = lines
            firstClique = True
            fastq+=1
        elif lines.startswith("@Clique2_"):
            sequences_paired1.write(buffer_lines.getvalue())
            buffer_lines = StringIO()
            name = lines
            firstClique = False
            fastq+=1
        elif firstClique:
            if fastq == 1:
                fastq+=1
                found = 0
                seq = lines
            elif fastq == 2:
                fastq+=1
            elif fastq == 3:
                fastq = 0
        else:
            if fastq == 1:
                found = 0
                fastq+=1
                seq = lines
            elif fastq == 2:
                fastq+=1
            elif fastq == 3:
                fastq = 0
        if lines.startswith("@Clique"):
            buffer_lines.write(lines.split("-+-")[0])
            buffer_lines.write("\n")
        else:
            buffer_lines.write(lines)

if firstClique:
    sequences_single.write(buffer_lines.getvalue())
    buffer_lines = StringIO()
else:
    sequences_paired2.write(buffer_lines.getvalue())
    buffer_lines = StringIO()

f = open("data_cliques_paired_R1.fastq",'w')
f.write(sequences_paired1.getvalue())
f.close()
f = open("data_cliques_paired_R2.fastq",'w')
f.write(sequences_paired2.getvalue())
f.close()
f = open("data_cliques_single.fastq",'w')
f.write(sequences_single.getvalue())
f.close()

count_buffer = {}
for genome in sys.argv[1:]:
    count_buffer[genome] = StringIO()

for singleDict in countlist:
    for genome in sys.argv[1:]:
        count = 0.0
        for read in singleDict:
            if genome in read:
                count += singleDict[read]/float(readnames[read])
                #print(singleDict[read],readnames[read])
        count_buffer[genome].write(str(count))
        count_buffer[genome].write("\n")

f = open("data_clique_to_reads.tsv", 'w')
f.write(clique_to_reads.getvalue())
f.close()
