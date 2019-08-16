#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/11 下午7:33
import sys
import os
import re
import argparse
from collections import defaultdict
from .Faidx import Faidx
from .helper import set_logging

SJMOTIF = set(["GT@AG", "CT@AC", "GC@AG", "CT@GC", "AT@AC", "GT@AT"])

LOG = set_logging("JunctionFetch")

def rev_strand(strand):
    """
    reverse the strand
    :param strand:
    :return:
    """
    if strand == "+":
        return "-"
    return "+"

class Junction:
    """
    A class to process the junction tab. All coordinary site don't need to convert to 0-based
    """

    def __init__(self, sj, genome, circ=False):

        allinfo = sj.split("_")
        self._name = sj
        self._info = sj[:-1].split("_")
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

    @staticmethod
    def checkmotif(motif, motifset):
        """
        :param motif:
        :param motifset:
        :return:
        """
        if motif in motifset:
            return True
        return False

    @staticmethod
    def junctionseq(faidx, chr, st, ed, strand, offset=0, circ=False):
        """

        :param faidx: genome Faidx
        :param chr: chromosome name
        :param st: start site
        :param ed: end site
        :param strand: plus or minus
        :param offset: for checking whether right site
        :return: upstream, downstream and combined sequence
        """
        assert isinstance(faidx, Faidx), "Faidx file must be a <class: Faidx>, not {}".format(faidx.__class__)
        if circ:
            seq1 = faidx.fetch(chr,
                               int(st), int(st) + 39,
                               rev_strand(strand))
            seq2 = faidx.fetch(chr,
                               int(ed) - 39, int(ed),
                               rev_strand(strand))
        else:
            seq1 = faidx.fetch(chr,
                               int(st) - 40 + offset,
                               int(st) - 1 + offset,
                               rev_strand(strand))

            seq2 = faidx.fetch(chr,
                               int(ed) - offset,
                               int(ed) + 39 - offset,
                               rev_strand(strand))

        if strand == "+":
            if circ:
                seq = seq2.realseq + seq1.realseq
            else:
                seq = seq1.realseq + seq2.realseq
        else:
            if circ:
                seq = seq1.realseq + seq2.realseq
            else:
                seq = seq2.realseq + seq1.realseq

        return seq1, seq2, seq

    @staticmethod
    def junctionmotif(faidx, chr, st, ed, strand, offset=0):
        """

        :param faidx:
        :param chr:
        :param st:
        :param ed:
        :param strand:
        :param offset:
        :return:
        """
        assert isinstance(faidx, Faidx), "Faidx file must be a <class: Faidx>, not {}".format(faidx.__class__)
        motif1 = faidx.fetch(chr, int(st) + offset, int(st) + 1 + offset, strand)
        motif2 = faidx.fetch(chr, int(ed) - 1 - offset, int(ed) - offset, strand)

        if strand == "+":
            motif = motif1.realseq + "@" + motif2.realseq
        else:
            motif = motif2.realseq + "@" + motif1.realseq

        return motif1, motif2, motif

    @staticmethod
    def circmotif(faidx, chr, st, ed, strand):
        """

        :param faidx:
        :param chr:
        :param st:
        :param ed:
        :param strand:
        :return:
        """
        assert isinstance(faidx, Faidx), "Faidx file must be a <class: Faidx>, not {}".format(faidx.__class__)
        motif1 = faidx.fetch(chr, int(st) - 1, int(st), strand)
        motif2 = faidx.fetch(chr, int(ed), int(ed) + 1, strand)

        if strand == "+":
            motif = motif2.realseq + "@" + motif1.realseq
        else:
            motif = motif1.realseq + "@" + motif2.realseq

        return motif1, motif2, motif

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
        u'''
        1.14 change the split object
        '''
        for j in self._info:

            u'''
            8.3 remove last string, that's the strand
            '''
            chr_, st, ed = re.split("-|:", j)

            if not self.circ:
                motif1, motif2, motif = Junction.junctionmotif(self.genome, chr_, st, ed, self.strand)
                seq1, seq2, seq = Junction.junctionseq(self.genome, chr_, st, ed, self.strand, circ=False)

                """
                We checked the splicing motif at the junction site. 
                if User gave a exon junction site, we would check the motif after correcting to intronic location.

                """

                if not Junction.checkmotif(motif, SJMOTIF):
                    motif1_, motif2_, motif_ = Junction.junctionmotif(self.genome, chr_, st, ed, self.strand, offset=1)
                    if Junction.checkmotif(motif_, SJMOTIF):
                        LOG.warn(
                            msg="{}\tIt seems the junction sites [{}] were on exon, not intronic. Already corrected.".format(
                                LOG.name, self._name))
                        motif1, motif2, motif = motif1_, motif2_, motif_
                        seq1, seq2, seq = Junction.junctionseq(self.genome, chr_, st, ed, self.strand, offset=1)

                res[j]["seq"] = seq
                res[j]["motif"] = motif
            else:
                motif1, motif2, motif = Junction.circmotif(self.genome, chr_, st, ed, self.strand)
                _, _, seq = Junction.junctionseq(self.genome, chr_, st, ed, self.strand, circ=True)

                res[j]["seq"] = seq
                res[j]["motif"] = motif

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
    LOG.info(msg="{}\tFetching the sequence from {} list.".format(LOG.name, junctionfile))
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
    configprefix = open(os.path.join(outdir, 'config.txt'), 'w')

    for j in junctionlist:
        junctionres = Junction(j, genome, circ=True)
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
            filelist.append(filename)
    configprefix.close()
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
