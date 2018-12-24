#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/11 下午7:33
import sys
import os
import re
import argparse
from collections import defaultdict
from Faidx import Faidx


class Junction:
    """
    A class to process the junction tab. All coordinary site don't need to convert to 0-based
    """

    def __init__(self, sj, genome, circ=False):

        allinfo = sj.split("_")
        self._name = sj
        self._info = sj.split("_")
        self._strand = allinfo[0][-1]
        self._chrom = allinfo[0].split(":")[0]
        self.circ = circ
        if isinstance(genome, Faidx):
            self._genome = genome
        else:
            self._genome = Faidx(genome)

    @property
    def chr(self):
        return self._chrom

    @property
    def strand(self):
        return self._strand

    @property
    def junctionname(self):
        return self._name

    @property
    def genome(self):
        return self._genome

    @property
    def juncmotif(self):
        """
        This is for return the canonical junction motif
        :return:
        """
        motif

    @property
    def junc_for_seq(self):
        """
        :param chr:
        :param junctionlist:
        :param strand:
        :return:
        """
        assert isinstance(self.genome, Faidx), "Faidx file must be a <class: Faidx>, not {}".format(faidx.__class__)
        res = defaultdict(lambda: defaultdict())
        for j in self._info:
            chr, st, ed = re.split("-|:", j[:-1])
            if not self.circ:
                seq1 = self.genome.fetch(chr, int(st) - 40, int(st) - 1, self.strand)
                motif1 = self.genome.fetch(chr, int(st), int(st) + 1, self.strand)
                seq2 = self.genome.fetch(chr, int(ed), int(ed) + 39, self.strand)
                motif2 = self.genome.fetch(chr, int(ed) - 1, int(ed), self.strand)
                if self.strand == "+":
                    seq = seq1.realseq + seq2.realseq
                else:
                    seq = seq2.realseq + seq1.realseq
                res[j]["seq"] = seq
                if self.strand == "+":
                    res[j]["motif"] = motif1.realseq + "@" + motif2.realseq
                else:
                    res[j]["motif"] = motif2.realseq + "@" + motif1.realseq
            else:
                seq1 = self.genome.fetch(chr, int(st), int(st) + 39, self.strand)
                motif1 = self.genome.fetch(chr, int(st) - 1, int(st), self.strand)
                seq2 = self.genome.fetch(chr, int(ed) - 40, int(ed) - 1, self.strand)
                motif2 = self.genome.fetch(chr, int(ed), int(ed) + 1, self.strand)
                if self.strand == "+":
                    seq = seq2.realseq + seq1.realseq
                else:
                    seq = seq2.realseq + seq1.realseq
                res[j]["seq"] = seq
                if self.strand == "+":
                    res[j]["motif"] = motif2.realseq + "@" + motif1.realseq
                else:
                    res[j]["motif"] = motif2.realseq + "@" + motif1.realseq
        return res


def _parse_args():
    userInput = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                        description="process the junction file and return the seq")

    userInput.add_argument('-fa', '--fasta', action='store', required=True,
                           type=str,
                           help='The genome file')
    userInput.add_argument('-ta', '--target', action='store', required=True,
                           type=str,
                           help='splicing junction file')
    userInput.add_argument('-od', '--outdir', action='store',
                           type=str, default="./",
                           help='output dir, defaul:./')
    userInput.add_argument('-circ', '--circRNA', action='store_true',
                           default=False,
                           help='circRNA, defaul:False')

    if len(sys.argv) == 1:
        userInput.print_help()
        sys.exit(1)
    return userInput.parse_args()


def checkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def targetfile(file):
    res = defaultdict(list)
    with open(file) as IN:
        for line in IN.readlines():
            line = line.strip().split("\t")
            res[line[0]] = line[1:]
            # res.append(line.strip().split("\t")[0])
    return res


def fetchjunc(genome, junctionfile, outdir):
    """
    :param args:
    :return:
    """

    filelist = []
    if outdir != "./":
        checkdir(outdir)
    genome = Faidx(genome)
    junctionlist = targetfile(junctionfile)
    configprefix = open(os.path.join(outdir, 'config.txt'), 'w')

    for j in junctionlist:
        junctionres = Junction(j, genome)
        index = 0
        for k, v in junctionres.junc_for_seq.items():
            filedir = os.path.join(outdir, junctionres.junctionname)
            newheader = ":".join([k, v["motif"]])
            filename = os.path.join(filedir, newheader + ".fasta")

            configprefix.write('\t'.join([newheader] + junctionlist[j][index].split(';')) + '\n')
            index += 1
            filelist.append(filename)
            checkdir(filedir)
            with open(filename, "w") as OUT:
                OUT.write(">" + newheader + "\n")
                OUT.write(v["seq"] + "\n")

    configprefix.close()
    return filelist


def fetchcirc(genome, junctionfile, outdir):
    filelist = []
    if outdir != "./":
        checkdir(outdir)
    genome = Faidx(genome)
    junctionlist = targetfile(junctionfile)
    for j in junctionlist:
        junctionres = Junction(j, genome, circ=True)
        for k, v in junctionres.junc_for_seq.items():
            filedir = os.path.join(outdir, junctionres.junctionname)
            newheader = ":".join([k, v["motif"]])
            filename = os.path.join(filedir, newheader + ".fasta")

            filelist.append(filename)
            checkdir(filedir)
            with open(filename, "w") as OUT:
                OUT.write(">" + newheader + "\n")
                OUT.write(v["seq"] + "\n")
            filelist.append(filename)
    return filelist


def main():
    # a = Junction('1:4837075-4839386+_1:4837075-4840955+',
    #              '/Volumes/bu15191450186/zr/singlecell/10X/refdata-cellranger-mm10-2.1.0/fasta/genome.fa')
    # for k, v in a.junc_for_seq().items():
    #     print(k, v.seq)
    args = _parse_args()
    if args.circRNA:
        fetchcirc(args.fasta, args.target, args.outdir)
    else:
        fetchjunc(args.fasta, args.target, args.outdir)


if __name__ == '__main__':
    main()
