#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/22 1:26 PM
__author__ = 'Zhou Ran'
import argparse
from GeneToTrans import generateinfo
from GenerateBlock import runSequenceCrawler
from AlignmentFilter import BlcokParser
from Bio.SeqUtils import MeltingTemp as mt

"""
这是第一个指令，利用基因ID来进行探针的设计的命令行整合
核心指令：
生成单个转录本序列和一个layer的配置文件信息。
1. python /Users/zhouran/opt/proj/2018-11-17-probe/probe/GeneToTrans.py  ../Mus_musculus.GRCm38.cdna.all.fa ../layer_specific_2.txt
2. python /Users/zhouran/opt/proj/2018-11-17-probe/probe/GenerateBlock.py -f ENSMUST00000020329.12.fasta -l 15 -L 22 -g 40 -G 60 -c 100 -C 300 -F 30 -O -o ENSMUST00000020329.12.block.fasta
3. python /Users/zhouran/opt/proj/2018-11-17-probe/probe/AlignmentFilter.py -x ../GRCm38/GRCm38_mod -s 390 -F 20 -T transfer.trans.txt -f ENSMUST00000020329.12.block.fasta.fastq -o a -v
"""


def transcript(args):
    fasta = args.fasta
    targetfile = args.target
    outputprefix = args.outprefix
    l = args.minLength
    L = args.maxLength
    gcPercent = args.min_GC
    GCPercent = args.max_GC
    tm = args.min_Tm
    TM = args.max_Tm
    X = args.prohibitedSeqs
    sal = args.salt
    form = args.formamide
    sp = args.Spacing
    concA = args.dnac1
    concB = args.dnac2
    headerVal = args.header
    bedVal = args.bed
    OverlapModeVal = args.OverlapMode
    verbocity = args.verbose
    reportVal = args.Report
    debugVal = args.Debug
    metaVal = args.Meta
    index = args.index

    if args.dnac1 >= args.dnac2:
        conc1 = args.dnac1
        conc2 = args.dnac2

    else:
        conc1 = args.dnac2
        conc2 = args.dnac1

    exec('nn_table = mt.%s' % args.nn_table, None, globals())
    falist = generateinfo(fasta, targetfile, outputprefix)
    for sub in falist:
        runSequenceCrawler(sub, l, L, gcPercent, GCPercent, nn_table, tm, TM,
                           X, sal, form, sp, conc1, conc2, headerVal, bedVal,
                           OverlapModeVal, verbocity, reportVal, debugVal, metaVal,
                           sub)
        BlcokParser('.'.join([sub, 'fastq']), index, '.'.join([outputprefix, 'layerinfo.txt']),
                    '.'.join([sub, 'result']), sal, form)


def arg():
    userInput = argparse.ArgumentParser(description= \
                                            '%s version %s. Requires a FASTA file as input. Currently, only '
                                            'single-entry FASTA files are supported.  Returns a .fastq file, which '
                                            'can be inputted into short read alignment programs. Optionally, a '
                                            '.bed file can be outputted instead if \'-b\' is flagged. Tm values '
                                            'are corrected for [Na+] and [formamide].' % ('aa', 'dd'))
    subparsers = userInput.add_subparsers(title='subcommands',
                                          description='valid subcommands',
                                          help='additional help',
                                          dest='subparser_name')

    parser_foo = subparsers.add_parser('trans')
    parser_foo.set_defaults(func=transcript)

    parser_foo.add_argument('-fa', '--fasta', action='store', required=True,
                            type=str,
                            help='The minimum allowed probe length; default is '
                                 '36')
    parser_foo.add_argument('-ta', '--target', action='store', required=True,
                            type=str,
                            help='The minimum allowed probe length; default is '
                                 '36')
    parser_foo.add_argument('-op', '--outprefix', action='store', required=True,
                            type=str,
                            help='The minimum allowed probe length; default is '
                                 '36')
    parser_foo.add_argument('-index', '--index', action='store', required=True,
                            type=str,
                            help='The minimum allowed probe length; default is '
                                 '36')
    parser_foo.add_argument('-l', '--minLength', action='store', default=36,
                            type=int,
                            help='The minimum allowed probe length; default is '
                                 '36')
    parser_foo.add_argument('-L', '--maxLength', action='store', default=41,
                            type=int,
                            help='The maximum allowed probe length, default is '
                                 '41')
    parser_foo.add_argument('-g', '--min_GC', action='store', default=20,
                            type=int,
                            help='The minimum allowed percent G + C, default is '
                                 '20')
    parser_foo.add_argument('-G', '--max_GC', action='store', default=80,
                            type=int,
                            help='The maximum allowed percent  G + C, default '
                                 'is 80')
    parser_foo.add_argument('-t', '--min_Tm', action='store', default=42,
                            type=int,
                            help='The minimum allowed Tm, default is 42')
    parser_foo.add_argument('-T', '--max_Tm', action='store', default=47,
                            type=int,
                            help='The maximum allowed Tm, default is 47')
    parser_foo.add_argument('-X', '--prohibitedSeqs', action='store',
                            default='AAAAA,TTTTT,CCCCC,GGGGG', type=str,
                            help='Prohibited sequence list (separated by commas '
                                 'with no spaces), default is '
                                 '\'AAAAA,TTTTT,CCCCC,GGGGG\'')
    parser_foo.add_argument('-s', '--salt', action='store', default=390,
                            type=int,
                            help='The mM Na+ concentration, default is 390')
    parser_foo.add_argument('-F', '--formamide', action='store', default=50,
                            type=float,
                            help='The percent formamide being used, default is '
                                 '50')
    parser_foo.add_argument('-S', '--Spacing', action='store', default=0,
                            type=int,
                            help='The minimum spacing between adjacent probes, '
                                 'default is 0 bases')
    parser_foo.add_argument('-c', '--dnac1', action='store', default=25,
                            type=float,
                            help='Concentration of higher concentration strand '
                                 '[nM] -typically the probe- to use for '
                                 'thermodynamic calculations. Default is 25')
    parser_foo.add_argument('-C', '--dnac2', action='store', default=25,
                            type=float,
                            help='Concentration of lower concentration strand '
                                 '[nM] -typically the target- to use for '
                                 'thermodynamic calculations. Default is 25')
    parser_foo.add_argument('-n', '--nn_table', action='store',
                            default='DNA_NN3',
                            type=str,
                            help='The nearest neighbor table of thermodynamic '
                                 'parameters to be used. See options in '
                                 'Bio.SeqUtils.MeltingTemp. Default is DNA_NN3')
    parser_foo.add_argument('-H', '--header', action='store', type=str,
                            help='Allows the use of a custom header in the '
                                 'format chr:start-stop. E.g. '
                                 '\'chr2:12500-13500\'')
    parser_foo.add_argument('-b', '--bed', action='store_true', default=False,
                            help='Output a .bed file of candidate probes '
                                 'instead of a .fastq file.')
    parser_foo.add_argument('-O', '--OverlapMode', action='store_true',
                            default=False,
                            help='Turn on Overlap Mode, which returns all '
                                 'possible candidate probes in a block of '
                                 'sequence including overlaps. Off by default. '
                                 'Note, if selecting this option, the '
                                 '-S/--Spacing value will be ignored')
    parser_foo.add_argument('-v', '--verbose', action='store_true',
                            default=False,
                            help='Turn on verbose mode to have probe mining'
                                 'progress print to Terminal. Off by default')
    parser_foo.add_argument('-R', '--Report', action='store_true', default=False,
                            help='Write a Report file detailing the results of '
                                 'each window of sequence considered by the '
                                 'script. The first set of lines give the '
                                 'occurrence of each possible failure mode for '
                                 'quick reference. Off by default. Note, '
                                 'selecting this option will slow the script '
                                 'considerably')
    parser_foo.add_argument('-D', '--Debug', action='store_true', default=False,
                            help='The same as -Report, but prints info to '
                                 'terminal instead of writing a log file. Off '
                                 'by default')
    parser_foo.add_argument('-M', '--Meta', action='store_true', default=False,
                            help='Write a text file containing meta information '
                                 'Off by default. Reports input file <tab> '
                                 'estimated runtime <tab> blockParse version '
                                 '<tab> candidates discovered <tab> span in kb '
                                 'covered by candidate probes <tab> candidate '
                                 'probes per kb')

    return userInput.parse_args()


def main():
    args = arg()
    args.func(args)


if __name__ == '__main__':
    main()
