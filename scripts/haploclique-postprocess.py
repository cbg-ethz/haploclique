#! /usr/bin/env python
import sys
from cStringIO import StringIO

firstClique = False
buffer_lines = StringIO()

fr1 = open("data_cliques_paired_R1.fastq",'w')
fr2 = open("data_cliques_paired_R2.fastq",'w')
fs = open("data_cliques_single.fastq",'w')
fc = open("data_clique_to_reads.tsv", 'w')

for lines in sys.stdin:
    if len(lines.rstrip()) > 0:
        if lines.startswith("@Clique1_"):
            if firstClique:
                fs.write(buffer_lines.getvalue())
                buffer_lines = StringIO()
            else:
                fr2.write(buffer_lines.getvalue())
                buffer_lines = StringIO()
            splits = lines[:-1].split("-+-")
            clique_name = splits[0][1:]
            fc.write("Clique_")
            fc.write(clique_name[8:])
            fc.write("\t")

            readSplit = splits[1][1:].split(",")
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
            start = 1
            for r in readCount:
                if r != "1" and r != "":
                    if (start == 0):
                        fc.write(",")
                    else:
                        start = 0
                    fc.write(str(readCount[r]))
                    fc.write("-")
                    fc.write(r)
            fc.write("\n")
            firstClique = True
        elif lines.startswith("@Clique2_"):
            fr1.write(buffer_lines.getvalue())
            buffer_lines = StringIO()
            firstClique = False

        if lines.startswith("@Clique"):
            buffer_lines.write(lines.split("-+-")[0])
            buffer_lines.write("\n")
        else:
            buffer_lines.write(lines)

if firstClique:
    fs.write(buffer_lines.getvalue())
    buffer_lines = StringIO()
else:
    fr2.write(buffer_lines.getvalue())
    buffer_lines = StringIO()

fr1.close()
fr2.close()
fs.close()
fc.close()
