#!/usr/bin/env python3

from subprocess import call
import random
import re
import argparse

class Sequence:

    def __init__(self, identifier, description):
        self.identifier = identifier
        self.description = description
        self.sequence = ''
        self.bases = ['A', 'G', 'C', 'T']

        self.base_freqs = {'A': 0, 'G': 0, 'C': 0, 'T': 0}

    def add_sequence(self, seq):
        old_len = len(self.sequence)
        self.sequence.extend(list(seq))

        length = len(self.sequence)

        new_len = len(seq)

        new_freq = {'A' : 0, 'G' : 0, 'C': 0, 'T': 0}

        for c in seq:
            new_freq[c] += 1

        for (base, freq) in new_freq.items():
            self.base_freqs[base] = (old_len * self.base_freqs[base] + new_len * freq) / length

    def sample_base(self):
        r = random.random()
        ct = 0

        for (base, freq) in new_freq.items():
            ct += freq
            if r < ct:
                return base

    def create_snp(self, pos=-1, base='N'):
        if pos < 0 or pos > len(self.sequence):
            pos = random.randint(0, len(self.ssequence)-1)

        if base == 'N':
            base = self.sample_base()

        self.sequence[pos] = base

    def create_insertion(self, mean, std_deviation, pos=-1):
        length = int(random.gauss(mean, std_deviation))

        if pos < 1 or pos > len(self.sequence):
            pos = random.randint(0, len(self.sequence)-1)

        insert = []

        for i in range(0, length):

            insert.append(self.sample_base())

        self.sequence[pos:pos] = insert

    def create_deletion(self, pos=-1, mean, std_deviation):
        length = int(random.gauss(mean, std_deviation))

        if pos < 1 or pos > len(self.sequence)-length:
            pos = random.randint(0, len(self.sequence)-length-1)

        self.sequence[pos:pos+length] = []

    def writeFASTA(filename):
        with open(filename, 'w') as f:
            f.write('>' + self.identifier + '|' + self.description + '\n')
            ct = 0
            for base in self.sequence:
                f.write(base)
                ct += 1
                if ct == 80:
                    ct = 0
                    f.write('\n')

class ParsingError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

def readFASTA(filename):
    with open(filename, 'r') as f:
        header = f.readline()

        m = re.search('>(?P<identifier>\w*)\|(?P<description>.*)', header)
        if not m:
            raise ParsingError('FASTA header was not valid')

        seq = Sequence(m.group('identifier'), m.group('description'))

        for line in f:
            seq.add_sequence(line)

    return seq

def main(argv):
    pass

if __name__ == '__main__':
    main(sys.argv)
