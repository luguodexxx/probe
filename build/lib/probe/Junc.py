#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/11 下午7:33

from Faidx import Faidx


# 1:21220132-21220132-21220534-21227158
# 1:21220132-21220235-21221843-21227158
# 1:21220132-21220132-21221746-21227158
# 1:21220132-21220534-21225441-21227158-21227227-21227930
# 1:21220132-21220845-21225441-21227158-21227227-21227930
# 1:21220132-21221746-21225441-21227158-21227227-21227930
# 1:21220132-21220235-21225441-21227158-21227227-21227930
# 1:21220132-21225249-21225441-21227158-21227227-21227930
# 1:21220132-21221746-21228042-21233557
# 1:21220132-21220845-21228042-21233557
# 1:21220132-21220845-21232106-21233557
# 1:21220132-21220534-21232106-21233557
# 1:21220132-21220132-21220845-21233557
# 1:21220132-21227158-21232106-21233557
# 1:21220132-21227158-21228042-21233557
# 1:21220132-21221746-21232106-21233557
# 1:21220132-21220534-21228042-21233557
# 1:21220132-21220132-21221746-21233557
# 1:21220220-21239900-21241728-21241728
class Junction:
    """
    A class to process the junction tab. All coordinary site don't need to convert to 0-based
    """

    def __init__(self, sj, genome):
        info = sj.split(':')
        jcsiteinfo = info[1].split('-')
        self._chrom = info[0]
        self._strand = info[-1]
        self._genome = Faidx(genome)
        if len(jcsiteinfo) == 4:
            self.s1, self.s2, self.s3, self.s4 = list(map(int, jcsiteinfo))
        else:
            self.s1, self.s2, self.s3, self.s4, self.s5, self.s6 = list(map(int, jcsiteinfo))

    @property
    def chr(self):
        return self._chrom

    @property
    def strand(self):
        return self._strand

    @property
    def genome(self):
        return self._genome

    @staticmethod
    def junc_for_seq(chr, junctionlist, strand, faidx):
        """

        :param chr:
        :param junctionlist:
        :param strand:
        :return:
        """
        assert isinstance(faidx, Faidx), "Faidx file must be a <class: Faidx>, not {}".format(faidx.__class__)
        res = []
        print(junctionlist)
        for j in junctionlist:
            seq1 = faidx.fetch(chr, j[0] - 40, j[0] - 1, strand)
            seq2 = faidx.fetch(chr, j[1], j[1] + 39, strand)
            if strand == "+":
                seq = seq1 + seq2
            else:
                seq = seq2 + seq1
            res.append(seq)
        return res


class SE(Junction):
    """
    A class to process SE events, three isoform specific sequence will return.
    """

    def __init__(self, sj, genome):
        Junction.__init__(self, sj, genome)

    @property
    def junc(self):
        return [[self.s1, self.s2], [self.s3, self.s4]], \
               [[self.s1, self.s4]]

    @property
    def seq(self):
        """
        The longest isoform was the first.
        :return:
        """
        iso1, iso2 = self.junc
        junc1 = Junction.junc_for_seq(self.chr, iso1, self.strand, self.genome)
        junc2 = Junction.junc_for_seq(self.chr, iso2, self.strand, self.genome)
        return [junc1, junc2]


class A3SS(Junction):
    """
    A class to process A3SS events, two isoform specific sequence will return.
    """

    def __init__(self, sj, genome):
        Junction.__init__(self, sj, genome)

    @property
    def junc(self):
        return [[self.s1, self.s3]], \
               [[self.s1, self.s4]]

    @property
    def seq(self):
        """
        The longest isoform was the first.
        :return:
        """
        iso1, iso2 = self.junc
        junc1 = Junction.junc_for_seq(self.chr, iso1, self.strand, self.genome)
        junc2 = Junction.junc_for_seq(self.chr, iso2, self.strand, self.genome)
        return [junc1, junc2]


class A5SS(Junction):
    """
    Same to A3SS
    """

    def __init__(self, sj, genome):
        Junction.__init__(self, sj, genome)

    @property
    def junc(self):
        return [[self.s1, self.s4]], \
               [[self.s2, self.s4]]

    @property
    def seq(self):
        """
        The longest isoform was the first.
        :return:
        """
        iso1, iso2 = self.junc
        junc1 = Junction.junc_for_seq(self.chr, iso1, self.strand, self.genome)
        junc2 = Junction.junc_for_seq(self.chr, iso2, self.strand, self.genome)
        return [junc1, junc2]


class MXE(Junction):
    """
    four isoform specific sequence will return
    """

    def __init__(self, sj, genome):
        Junction.__init__(self, sj, genome)

    @property
    def junc(self):
        return [[self.s1, self.s2], [self.s3, self.s6]], \
               [[self.s1, self.s4], [self.s5, self.s6]]

    @property
    def seq(self):
        """
        The longest isoform was the first.
        :return:
        """
        iso1, iso2 = self.junc
        junc1 = Junction.junc_for_seq(self.chr, iso1, self.strand, self.genome)
        junc2 = Junction.junc_for_seq(self.chr, iso2, self.strand, self.genome)
        return [junc1, junc2]


class SingleJ:
    def __init__(self, sj, genome):
        info = sj.split(':')
        self._chr = info[0]
        self.s1, self.s2 = info[1].split('-')
        self.strand = info[2]
        self.genome = Faidx(genome)

    @property
    def junc(self):
        return [self.s1, self.s2]

    @property
    def seq(self):
        junc = Junction.junc_for_seq(self._chr, self.junc, self.strand, self.genome)
        return junc


def main():
    a = A5SS('1:21220132-21221746-21228042-21233557',
             '/mnt/raid61/Microwell/mm10/fasta/genome.fa')
    # genomefa = Faidx('/Users/zhouran/opt/proj/2018-11-17-probe/rebuild/Probes_evaluation/testFaidx.fa')
    print(a.junc)
    b = a.seq
    print(b)


if __name__ == '__main__':
    main()
