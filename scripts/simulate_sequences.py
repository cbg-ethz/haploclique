#!/usr/bin/env python3

from subprocess import call
import random
import re
import argparse
import tarfile
import tempfile

class Sequence:

    def __init__(self, identifier, description):
        self.identifier = identifier
        self.description = description
        self.sequence = ''
        self.bases = ['A', 'G', 'C', 'T']

        self.base_freqs = {'A': 0, 'G': 0, 'C': 0, 'T': 0}

    def replicate(self):
        new_seq = Sequence(self.identifier, self.description)
        new_seq.add_sequence(''.join(self.sequence))

        return new_seq

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

    def writeFASTA(fileobj):
        fileobj.write('>' + self.identifier + '|' + self.description + '\n')
        ct = 0
        for base in self.sequence:
            fileobj.write(base)
            ct += 1
            if ct == 80:
                ct = 0
                fileobj.write('\n')

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
    parser = argparse.ArgumentParser(description='Generate haplotypes from a reference genome and create simulated illumina reads with SimSeq')
    parser.add_argument('-i', '--input', help='the reference genome in FASTA format (REQUIRED)', metavar='reference.fasta')
    parser.add_argument('-o', '--output', help='destination of the generated data tar.gz archive. It contains the reference sequence, the generated haplotypes and the simulated sam file.', metavar='output.tar.gz')
    parser.add_argument('--snp', help='generate new haplotypes through substituting random nucleotides.')
    parser.add_argument('--ins', help='generate new haplotypes through inserting random sequences.')
    parser.add_argument('--del', help='generate new haplotypes through deleting random sequences.')
    parser.add_argument('-m', '--mean', type=int, help='mean value for length of indels or number of nucleotide substitutions')
    parser.add_argument('-s', '--sigma', type=int, help='standard deviation for length of indels or number of nucleotide substitions')
    parser.add_argument('-n', '--number', type=int, help='number of created haplotypes')
    parser.add_argument('--seed', type=int, help='seed for the Random Number Generator')

    args = parser.parse_args(argv)

    ref_genome = readFASTA(args.input)

    m = re.match('(?P<name>.*)(\.fasta$)', args.input)
    filename = m.group('name')

    if not args.output:
        args.output = filename + '_sht.tar.gz'

    tar = tarfile.open(args.output, 'w:gz')
    tar.add(args.input, arcname='ref_' + args.input)

    n = args.n

    haplofiles = []

    for i in range(0, n):
        seq = ref_genome.replicate()

        mean = args.mean
        sigma = args.sigma

        if mean < 1:
            mean *= len(ref_genome.sequence)
            sigma *= len(ref_genome.sigma)

        if args.snp:
            for i in range(0, random.gauss(mean, sigma))
                seq.create_snp()

        if args.ins:
            seq.create_insertion(mean, sigma)

        if args.del:
            seq.create_deletion(mean, sigma)

        f = tempfile.NamedTemporaryFile()
        seq.writeFASTA(f)
        info = tar.gettarinfo(arcname = 'ht' + i + '_' + filename + '.fasta',fileobj=f)
        tar.addfile(info, f)

        haplofiles.append(f)

    for hf in haplofiles:

        call(['../bin/SimSeq.jar', #TODO 
        hf.close()

    tar.close()

if __name__ == '__main__':
    main(sys.argv)
