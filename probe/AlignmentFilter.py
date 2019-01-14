#!/usr/bin/env python

import sys
import os
import glob
from collections import defaultdict
from itertools import chain
import random
import subprocess
import optparse
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import multiprocessing
from .helper import set_logging

SJMOTIF = set(["GT@AG", "CT@AC", "GC@AG", "CT@GC", "AT@AC", "GT@AT"])
LOG = set_logging("AlignmentFilter")


class Bowtie2Error(Exception):
    pass


# def mfold(falist, ct, na_conc, type="DNA"):
def mfold(falist, ct, na_conc):
    """
    calling the mfold and check the second structure
    :param fa:
    :param type:
    :param NA_CONC:
    :param Tm:
    :return:
    """
    # falist, ct, na_conc = args
    faprefix, left, right = falist[0].split(';')

    with open(faprefix + '.fa', 'w') as Fa:
        Fa.write(">{}\n{}".format(falist[0], falist[5]))

    mfoldcomm = "mfold_mod SEQ=\'{}\' NA={} NA_CONC={} T={}".format(faprefix + '.fa', "DNA", na_conc, int(ct))

    proc = subprocess.Popen(
        mfoldcomm,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )

    result, err = proc.communicate()

    try:
        with open(faprefix + ".fa.ct") as IN:
            res = []
            for line in IN.readlines():
                res.extend([line.strip().split()])
            leftcheck = set(map(lambda x: x[4], res[1:int(right)]))
            rightcheck = set(map(lambda x: x[4], res[-int(left):]))
    except:
        return False

    for f in glob.glob("{}.fa*".format(faprefix)):
        # print(f)
        os.unlink(f)

    if len(leftcheck) == 1 and len(rightcheck) == 1:
        return falist
    else:
        return False


def wrapperprocess(args):
    """

    :param args:
    :return:
    """
    return mfold(*args)


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
        if not self.mapped:
            self.geneid, self.genesymbol, self.location = ["nogeneid", "nogenesymbol", "0:0:0:1"]
        else:
            self.geneid, self.genesymbol, self.location = self.geneinfo.split("|")

    @property
    def proRC(self):
        return str(Seq(self.seq, IUPAC.unambiguous_dna).reverse_complement())

    @property
    def PLP(self):
        left, right = self.f_chr.split(';')[1:]

        return (self.proRC[:int(left)], self.proRC[int(left):])

    @property
    def mapped(self):
        if int(self.flag) & 4:
            return False
        return True

    @property
    def internal_start(self):
        return self.location.split(":")[1]

    @property
    def internal_end(self):
        return self.location.split(":")[2]

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
        Check the reads were mapped on minus[True] or plus[False]
        """
        if int(self.flag) & 16 == 16:
            return True
        else:
            return False

    def __str__(self):
        return self.info


def locationfilter(fake_chrname):
    """
    Filter the probe location. default length: 80nt
    :param junctionname: @chrom:2-39;19;19
    :return:
    """
    interval, uplen, downlen = fake_chrname.split(';')
    upstream, downstream = interval.split(':')[1].split('-')

    if int(upstream) <= 36 and int(downstream) >= 45:
        return True
    else:
        return False


class JuncParser():
    """
    JunctParser
    """

    def __init__(self,
                 fa,
                 index,
                 targetfile,
                 outfile,
                 sal,
                 formamide,
                 probelength,
                 hytemp,
                 thread,
                 mfold_=False,
                 verbose=False):
        self.fa = fa
        self._prefix = os.path.splitext(os.path.split(self.fa)[1])[0]
        self.st, self.ed = self._prefix.split(":")[1][:-1].split('-')
        self.index = index
        self._targetfile = targetfile
        self.TARGET = BlockParser.processtarget(self._targetfile)
        self.sal = sal
        self.verbose = verbose
        self.formamide = formamide
        self._probelength = probelength
        self.outfile = outfile
        self.hytemp = hytemp
        self.correcttemp = 0.65 * self.formamide + self.hytemp
        self.mfold = mfold_
        self.samresult = BlockParser.processAlign(self.index, self.fa, self.sal, self.formamide)
        self.filter = self.__filter()

        pool = multiprocessing.Pool(processes=thread)
        results = []
        if self.verbose:
            LOG.info(msg="{}\tWriting the bowtie2 results to {}".format(LOG.name, self.outfile + '.sam'))
            with open(os.path.splitext(self.outfile)[0] + '.sam', 'w') as tmpSAM:
                for read in self.samresult:
                    tmpSAM.write(str(read) + '\n')

        LOG.info(msg="{}\tWriting the results to {}".format(LOG.name, self.outfile))

        with open(self.outfile, 'w') as OUT:
            OUT.write('\t'.join(
                ['FakeChrom', 'motif', 'canonical', 'left', 'right', 'afterRC', 'beforeRC',
                 "PLPsequence", 'Tm', 'isoform_nums', 'isoforms']) + '\n')
            for read in self.filter:
                if self.mfold:
                    args = (read, self.correcttemp, self.sal / 1000)
                    results.append(pool.apply_async(wrapperprocess, args=(args,)))
                else:
                    OUT.write('\t'.join(read) + '\n')

            if self.mfold:
                pool.close()
                pool.join()
                os.system('rm *fa* -rf')
                for res in results:
                    res = res.get()
                    if res:
                        OUT.write('\t'.join(res) + '\n')

    @property
    def motif(self):
        return self._prefix.split(':')[-1]

    def __filter(self):
        """
        :return:
        """

        samdict = defaultdict(lambda: defaultdict(list))
        probeseqinfo = self.TARGET[self._prefix][:3]

        result = []
        for line in self.samresult:
            if not line.mapped:
                chrom, start, stop, seq, Tm, revseq = line.f_chr, line.abs_start, line.abs_end, line.seq, line.Tm, \
                                                      line.proRC
                left, right = line.PLP
                plpseq = generateprobe(left, right, self._probelength, probeseqinfo, self._prefix)

                result.append((chrom, self.motif,
                               "yes" if self.motif in SJMOTIF else "no", left, right, revseq, seq, plpseq, Tm, '1',
                               line.geneinfo))
                continue
            samdict[line.f_chr][line.location].append(line)

        for fakechrom, internalinfo in samdict.items():
            if locationfilter(fakechrom):
                for location, samline in internalinfo.items():
                    if samline.internal_start <= self.st and samline.internal_end >= self.ed:
                        result.append(line)
                    else:
                        if not line.checkFlag():
                            chrom, start, stop, seq, Tm, revseq = line.f_chr, line.abs_start, line.abs_end, line.seq, line.Tm, \
                                                                  line.proRC
                            left, right = line.PLP
                            plpseq = generateprobe(left, right, self._probelength, probeseqinfo, self._prefix)

                            result.append(
                                (chrom, self.motif, "yes" if self.motif in SJMOTIF else "no", left, right, revseq, seq,
                                 plpseq, Tm, '1', line.geneinfo))
                        else:
                            continue
        return result


class BlockParser():
    """
    BlockParser, execute the bowtie2 command and filter every mapped information on the air
    """

    def __init__(self, fa, index, tagetfile, outfile, sal, formamide, probelength, hytemp, thread, mfold_=False,
                 verbose=False):
        self.fa = fa
        self._prefix = os.path.splitext(os.path.split(self.fa)[1])[0]
        self.index = index
        self._targetfile = tagetfile

        self.TARGET = BlockParser.processtarget(self._targetfile)
        self.outfile = outfile
        self.sal = sal
        self.verbose = verbose
        self.formamide = formamide
        self._probelength = probelength

        self.samresult = BlockParser.processAlign(self.index, self.fa, self.sal, self.formamide)

        self.hytemp = hytemp
        self.correcttemp = 0.65 * self.formamide + self.hytemp
        self.mfold = mfold_
        self.filter = self.__filter()

        pool = multiprocessing.Pool(processes=thread)
        results = []
        if self.verbose:

            LOG.info(msg="{}\tWriting the bowtie2 results to {}".format(LOG.name, self.outfile + '.sam'))
            with open(self.outfile + '.sam', 'w') as tmpSAM:
                for read in self.samresult:
                    tmpSAM.write(str(read) + '\n')

        LOG.info(msg="{}\tWriting the results to {}".format(LOG.name, self.outfile))
        with open(self.outfile, 'w') as OUT:
            OUT.write('\t'.join(
                ['FakeChrom', 'left', 'right', 'afterRC', 'beforeRC',
                 "PLPsequence", 'Tm', 'isoform_nums', 'isoforms', 'symbolId']) + '\n')
            for read in self.filter:
                if self.mfold:
                    args = (read, self.correcttemp, self.sal / 1000)
                    results.append(pool.apply_async(wrapperprocess, args=(args,)))
                else:
                    OUT.write('\t'.join(read) + '\n')

            if self.mfold:
                pool.close()
                pool.join()
                os.system('rm *fa* -rf')
                for res in results:
                    res = res.get()
                    if res:
                        OUT.write('\t'.join(res) + '\n')

    @staticmethod
    def processtarget(targetfile):
        targetres = defaultdict(list)
        with open(targetfile) as IN:
            for line in IN.readlines():
                line = line.strip().split('\t')
                targetres[line[0]] = line[1:]
        return targetres

    @staticmethod
    def process_revline(samlist):
        """
        if a candidate only aligned to the host gene, we only choose the reverse alignment
        :param samlist:
        :return:
        """
        transcriptid = set()
        for rline in samlist:
            if rline.checkFlag():
                transcriptid.add('|'.join([rline.geneid, rline.genesymbol]))
            else:
                continue
        return transcriptid

    @staticmethod
    def process_revline_multihostgene(samlist):
        """
        if a candidate aligned to multiple genes and including the host gene,
        not reverse alignment on other genes also accepted, accepted[True] not[False]
        :param samlist:
        :return:
        """
        results = set()
        for rline in samlist:
            if rline.checkFlag():
                return False, results
            results.add('|'.join([rline.geneid, rline.genesymbol]))

        return True, results

    @staticmethod
    def processAlign(index, fa, sal, formamide):
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
            '-x', index,
            '-U', fa,
            params
        ]

        LOG.info(msg='{}\tBowtie2 command: {}'.format(LOG.name, ' '.join(bowtie2comm)))

        proc = subprocess.Popen(
            bowtie2comm,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        result, err = proc.communicate()
        if not result:
            err = err.decode()
            if "ERR" in err:
                raise Bowtie2Error("{}".format(err))

        parsing = [
            SAM(line, sal, formamide)
            for line in result.decode('utf-8').splitlines()
        ]

        return parsing

    def __filter(self):
        """
        Process the mapped information after mapped
        """

        result = []
        probeseqinfo = self.TARGET[self._prefix][:3]
        genesymbol = self.TARGET[self._prefix][4]
        additional = self.TARGET[self._prefix][3]

        samdict = defaultdict(lambda: defaultdict(list))
        for line in self.samresult:
            samdict[line.f_chr][line.genesymbol].append(line)

        for fakeChr, genesymbol_ref in samdict.items():
            referkeys = list(genesymbol_ref.keys())
            if len(referkeys) == 1:
                if genesymbol == referkeys[0]:
                    line_ = list(chain(*genesymbol_ref.values()))
                    transcriptid = BlockParser.process_revline(line_)
                    chrom, seq, Tm, revseq = line_[0].f_chr, line_[0].seq, line_[0].Tm, \
                                             line_[0].proRC
                    left, right = line_[0].PLP
                    plpseq = generateprobe(left, right, self._probelength, probeseqinfo, self._prefix)
                    result.append(
                        (chrom, left, right, revseq, seq, plpseq, Tm, str(len(transcriptid)), ','.join(transcriptid),
                         additional)
                    )
                else:
                    continue
            else:
                if genesymbol in set(referkeys):
                    linelist = list(chain(*[v for k, v in genesymbol_ref.items() if k != genesymbol]))
                    passstatus, reslist = BlockParser.process_revline_multihostgene(linelist)
                    if passstatus:
                        line_ = genesymbol_ref[genesymbol]
                        transcriptid = list(BlockParser.process_revline(line_)) + list(reslist)
                        chrom, seq, Tm, revseq = line_[0].f_chr, line_[0].seq, line_[0].Tm, \
                                                 line_[0].proRC
                        left, right = line_[0].PLP
                        plpseq = generateprobe(left, right, self._probelength, probeseqinfo, self._prefix)

                        result.append(
                            (chrom, left, right, revseq, seq, plpseq, Tm, str(len(transcriptid)),
                             ','.join(transcriptid), additional)
                        )
                else:
                    continue
        return result


def generateprobe(left, right, probelength, configinfo, hostname, gccontent=0.5):
    """

    :param samline:
    :param configinfo:
    :return:
    """
    firbc, secbc, thirdbc = configinfo
    alignedlength = len(left + right)
    retain = probelength - (alignedlength + len(firbc) + len(secbc) + len(thirdbc))
    if retain < 0:
        LOG.warn(
            "probelength was too short. all: {}, aligned: {}, fir,sec,third {} {} {}. HostGene: {}".format(probelength,
                                                                                                           alignedlength,
                                                                                                           len(firbc),
                                                                                                           len(secbc),
                                                                                                           len(thirdbc),
                                                                                                           hostname))
        return " "
    DNA = ["A", "T",
           "C", "G"]

    leftrandom = int(retain / 2)
    rightrandom = retain - leftrandom
    leftgc = [0 for i in range(int(leftrandom * gccontent))] + \
             [2 for i in range(int(leftrandom - int(leftrandom * gccontent)))]

    rightgc = [0 for i in range(int(rightrandom * gccontent))] + \
              [2 for i in range(int(rightrandom - int(rightrandom * gccontent)))]

    random.shuffle(leftgc)
    random.shuffle(rightgc)

    leftseq = ''.join([DNA[random.choice([0, 1]) + i] for i in leftgc])
    rightseq = ''.join([DNA[random.choice([0, 1]) + i] for i in rightgc])

    # leftseq = ''.join([random.choice(['A', 'T', 'C', 'G']) for i in range(leftrandom)])
    # rightseq = ''.join([random.choice(['A', 'T', 'C', 'G']) for i in range(rightrandom)])

    return "".join([right, leftseq, firbc, secbc, thirdbc, rightseq, left])


def _parse_args():
    """
    Parse the command line for options
    """
    usage = """AlignmentFilter.py -x [bowtie2index] -f [FqFile] -s <390> -F <20> -T [TargetPoolFile] -o <output>
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
        '-p', '--probelength', action='store', default=70, type=int, dest='probelength',
        help='probelength.'
    )
    group.add_option(
        '-h', '--hp', action='store', default=47, type=float, dest='hytemp',
        help="hybridize temp"
    )

    parser.add_option_group(group)
    options, args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return options


def main():
    options = _parse_args()
    samresult = BlockParser(
        options.file,
        options.index,
        options.targets,
        options.output,
        options.salt,
        options.formamide,
        options.probelength,
        options.hytemp,
        options.verbose
    )
    # samresult.__filter()


if __name__ == '__main__':
    main()
