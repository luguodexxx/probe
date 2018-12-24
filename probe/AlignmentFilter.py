#!/usr/bin/env python
# --------------------------------------------------------------------------
#
# (c) 2018 ChenLin's Lab
# 
# shengwuzhiliaoguojiazhongdianshiyanshi
# Sichuan University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# --------------------------------------------------------------------------


import sys
import os
from collections import defaultdict
import random
import subprocess
import optparse
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class IndexError(Exception):
    pass


class SAM():
    """
    Process the every mapped information
    """

    def __init__(self, line, sal, form):
        # f_ means fake, for exampel, f_chr -> fake chromsome id
        self.info = line  # contain all information
        self.f_chr, self.flag, self.geneinfo, self.abs_start, self.abs_end, self.cigar = line.split('\t')[:6]
        self.seq = line.split('\t')[9]
        self.sal = int(sal)
        self.form = form
        self.Tm = self.__probeTm()

    @property
    def geneid(self):
        return self.geneinfo.split("_")[0]

    @property
    def proRC(self):
        return str(Seq(self.seq, IUPAC.unambiguous_dna).reverse_complement())

    @property
    def PLP(self):
        # print(self.f_chr)
        left, right = self.f_chr.split(';')[1:]

        return (self.proRC[:int(left)], self.proRC[int(left):])

    def __probeTm(self):
        """
        Calculates the melting temperature of a given sequence under the
        specified salt and formamide conditions.
        """
        tmval = float(mt.Tm_NN(self.seq, Na=self.sal))
        Tm = ('%0.2f' % mt.chem_correction(tmval, fmd=int(self.form)))
        return Tm

    def checkFlag(self):
        """
        Check the reads were mapped on plus[True] or minus[False]
        """
        if int(self.flag) & 16 == 16:
            return True
        else:
            return False

    def __str__(self):
        return self.info


class BlcokParser():
    """
    BlcokParser, execute the bowtie2 command and filter every mapped information on the air
    """

    def __init__(self, fa, INDEX, TARGETfile, OUTfile, sal, formamide, probelength, verbose=False):
        self.fa = fa
        self._prefix = os.path.splitext(self.fa)[0]
        self.index = INDEX
        self._targetfile = TARGETfile
        # self.TARGET = set([i.strip().split('\t')[0] for i in open(TARGETfile).readlines()])
        self.TARGET = self.__processtarget()
        self.OUTfile = OUTfile
        self.sal = sal
        self.verbose = verbose
        self.formamide = formamide
        self._probelength = probelength

        self.samresult = self.__processAlign()
        self.filter = self.__filter()

        if self.verbose:
            with open('tmp.sam', 'w') as tmpSAM:
                for read in self.samresult:
                    tmpSAM.write(str(read) + '\n')

        with open(self.OUTfile, 'w') as OUT:
            OUT.write('\t'.join(
                ['FakeChrom', 'left', 'right', 'afterRC', 'beforeRC', "PLPsequence", 'Tm', 'isoform_nums',
                 'isoforms']) + '\n')
            for read in self.filter:
                OUT.write('\t'.join(read) + '\n')

    def __processtarget(self):
        targetres = defaultdict(list)
        with open(self._targetfile) as IN:
            for line in IN.readlines():
                line = line.strip().split('\t')
                targetres[line[0]] = line[1:]
        return targetres

    def __processAlign(self):
        """
        Process the bowtie command on the air
        """
        checkcomm = subprocess.call(
            ['which', 'bowtie2'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        if checkcomm != 0:
            raise FileNotFoundError('bowtie2 were not found! please install bowtie2')

        bowtie2 = 'bowtie2'
        params = '--no-hd -t -k 10 --local -D 20 -R 3 -N 1 -L 20 -i C,4 --score-min G,1,4 --end-to-end'
        bowtie2comm = [
            bowtie2,
            '-x', self.index,
            '-U', self.fa,
            params
        ]

        print('Bowtie2 command:\n {}'.format(' '.join(bowtie2comm)))
        proc = subprocess.Popen(
            bowtie2comm,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        result, err = proc.communicate()
        if not result:
            raise IndexError("Could not locate a Bowtie index corresponding to basename \"{}\"".format(self.index))

        parsing = [
            SAM(line, self.sal, self.formamide)
            for line in result.decode('utf-8').splitlines()
        ]

        return parsing

    def __filter(self):
        """
        Process the mapped information after mapped
        """

        result = []
        probeseqinfo = self.TARGET[self._prefix][:3]
        DIC = defaultdict(lambda: defaultdict(list))
        for line in self.samresult:
            DIC[line.f_chr][line.geneid].append(line)

        for fakeChr, GeneId in DIC.items():
            alignment_target = GeneId.keys()
            if any(sub_gene != self._prefix for sub_gene in alignment_target):
                A = []  # in TARGET
                B = []  # out of TARGET
                for gene in alignment_target:
                    if gene == self._prefix:
                        for line in GeneId[gene]:
                            if int(line.flag) & 16 == 16:
                                A.append(line)
                    else:
                        for line in GeneId[gene]:
                            B.append(line)
                if any(i.checkFlag() for i in B) and len(A) != 0:
                    isoform = [saminfo.geneid for saminfo in A]
                    chrom, start, stop, seq, Tm, revseq = A[0].f_chr, A[0].abs_start, A[0].abs_end, A[0].seq, A[0].Tm, \
                                                          A[0].proRC
                    left, right = A[0].PLP
                    plpseq = generateprobe(left, right, self._probelength, probeseqinfo)
                    result.append((chrom, left, right, revseq, seq, plpseq, Tm, str(len(isoform)), ','.join(isoform)))
            else:
                res = []
                for gene, saminfo in GeneId.items():
                    for i in saminfo:
                        if i.checkFlag():
                            res.append(i.geneid)
                            chrom, start, stop, seq, Tm, revseq = i.f_chr, i.abs_start, i.abs_end, i.seq, i.Tm, i.proRC
                            left, right = i.PLP
                plpseq = generateprobe(left, right, self._probelength, probeseqinfo)
                result.append((chrom, left, right, revseq, seq, plpseq, Tm, str(len(res)), ','.join(res)))
        return result


def generateprobe(left, right, probelength, configinfo):
    """

    :param samline:
    :param configinfo:
    :return:
    """
    firbc, secbc, thirdbc = configinfo
    alignedlength = len(left + right)
    retain = probelength - (alignedlength + len(firbc) + len(secbc) + len(thirdbc))
    assert retain >= 0, "probelength was too short. all: {}, aligned: {}, fir,sec,third {} {} {}".format(probelength,
                                                                                                         alignedlength,
                                                                                                         len(firbc),
                                                                                                         len(secbc),
                                                                                                         len(thirdbc))
    leftrandom = int(retain / 2)
    rightrandom = retain - leftrandom
    leftseq = ''.join([random.choice(['A', 'T', 'C', 'G']) for i in range(leftrandom)])
    rightseq = ''.join([random.choice(['A', 'T', 'C', 'G']) for i in range(rightrandom)])

    return "".join([right, leftseq, firbc, secbc, thirdbc, rightseq, left])


def _parse_args():
    """
    Parse the command line for options
    """
    usage = """callpeaks.py -x [bowtie2index] -f [FqFile] -s <390> -F <20> -T [TargetPoolFile] -o <output>
            """
    fmt = optparse.IndentedHelpFormatter(max_help_position=50, width=100)
    parser = optparse.OptionParser(usage=usage, formatter=fmt)
    group = optparse.OptionGroup(parser, 'Query arguments',
                                 'These options define search query arguments and parameters.                                    '
                                 'bowtie2 command has been the optimal params.')

    group.add_option(
        '-f', '--file', action='store', dest='file', type='str',
        help='The fastq file to be processed'
    )

    group.add_option(
        '-s', '--salt', action='store', default=390, type=int, dest='salt',
        help='The mM Na+ concentration to be used for Tm '
             'calculation, default is 390'
    )

    group.add_option(
        '-T', '--targets', action='store', dest='targets',
        help='Excluding these reads were mapped on non-target genes'
    )

    group.add_option(
        '-F', '--formamide', action='store', default=50, type=float, dest='formamide',
        help='The percent formamide to be used for Tm '
             'calculation, default is 50'
    )

    group.add_option(
        '-o', '--output', action='store', default=None, type=str, dest='output',
        help='Specify the stem of the output filename'
    )

    group.add_option(
        '-x', '--index', action='store', default=None, type=str, dest='index',
        help='The genome index dirs which generated by bowtie2'
    )

    group.add_option(
        '-v', '--verbose', action='store_true', default=False, dest='verbose',
        help='if True, the bowtie2 alignment results will stroe in the tmp.sam file.default:False'
    )
    group.add_option(
        '-p', '--probelength', action='store', default=None, type=int, dest='probelength',
        help='probelength.'
    )

    parser.add_option_group(group)
    options, args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return options


def main():
    options = _parse_args()
    samresult = BlcokParser(
        options.file,
        options.index,
        options.targets,
        options.output,
        options.salt,
        options.formamide,
        options.probelength,
        options.verbose
    )


if __name__ == '__main__':
    main()
