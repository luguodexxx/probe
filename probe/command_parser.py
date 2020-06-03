#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/22 1:26 PM
__author__ = 'Zhou Ran'

import os
import sys
import argparse
import textwrap
from Bio.SeqUtils import MeltingTemp as mt

from .GeneToTrans import generateinfo
from .GenerateBlock import runSequenceCrawler
from .AlignmentFilter import BlockParser, JuncParser
from .Junc import fetchjunc, fetchcirc
from .version import __version__


def transcript(args):
    """
    too many parameters!!!
    :param args:
    :return:
    """
    fasta = args.fastaC
    # tmp = open(args.config).read()
    # print(tmp)
    exec(open(args.config).read(), None, globals())
    exec('nn_table = mt.%s' % nn_table, None, globals())
    falist = generateinfo(fasta, targetfile, outputprefix)
    for sub in falist:
        subprefix = os.path.splitext(sub)[0]
        runSequenceCrawler(sub, l, L, gcPercent, GCPercent, nn_table, tm, TM, X, sal, form, sp, conc1, conc2,
                           OverlapModeVal, subprefix, entropy, vargibbs, par, saltscheme, ct)

        BlockParser('.'.join([subprefix, 'fastq']), index, '.'.join([outputprefix, 'layerinfo.txt']),
                    '.'.join([subprefix, 'result']), sal, form, probelength, hytemp, thread, detG, cDNA, mfold_=mfold,
                    verbose=verbocity)


def junction(args):

    fastaG = args.fastaG
    exec(open(args.config).read(), None, globals())
    exec('nn_table = mt.%s' % nn_table, None, globals())
    falist = fetchjunc(fastaG, targetfile, outputprefix)
    for sub in falist:
        subprefix = os.path.splitext(sub)[0]
        runSequenceCrawler(sub, l, L, gcPercent, GCPercent, nn_table, tm, TM, X, sal, form, sp, conc1, conc2,
                           OverlapModeVal, subprefix, entropy,vargibbs, par, saltscheme, ct)

        JuncParser('.'.join([subprefix, 'fastq']), index, os.path.join(outputprefix, 'config.txt'),
                   '.'.join([subprefix, 'result']), sal, form, probelength, hytemp, thread, detG, cDNA, mfold_=mfold,
                   verbose=verbocity)


def circ(args):
    fastaG = args.fastaG
    exec(open(args.config).read(), None, globals())
    exec('nn_table = mt.%s' % nn_table, None, globals())

    falist = fetchcirc(fastaG, targetfile, outputprefix)
    for sub in falist:
        subprefix = os.path.splitext(sub)[0]
        runSequenceCrawler(sub, l, L, gcPercent, GCPercent, nn_table, tm, TM, X, sal, form, sp, conc1, conc2,
                           OverlapModeVal, subprefix, entropy, vargibbs, par, saltscheme, ct)

        JuncParser('.'.join([subprefix, 'fastq']), index, os.path.join(outputprefix, 'config.txt'),
                   '.'.join([subprefix, 'result']), sal, form, probelength, hytemp, thread, detG, cDNA, mfold_=mfold,
                   verbose=verbocity)


def arg():
    probedesign = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False,
                                          description=textwrap.dedent("""
__________              ___.          ________                .__               
\______   \_______  ____\_ |__   ____ \______ \   ____   _____|__| ____   ____  
 |     ___/\_  __ \/  _ \| __ \_/ __ \ |    |  \_/ __ \ /  ___/  |/ ___\ /    \ 
 |    |     |  | \(  <_> ) \_\ \  ___/ |    `   \  ___/ \___ \|  / /_/  >   |  \

 |____|     |__|   \____/|___  /\___  >_______  /\___  >____  >__\___  /|___|  /
                             \/     \/        \/     \/     \/  /_____/      \/                              """))

    probedesign.add_argument('-config', '--config', action='store',
                             type=str, required=True,
                             help='A config file')

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""
__________              ___.          ________                .__               
\______   \_______  ____\_ |__   ____ \______ \   ____   _____|__| ____   ____  
 |     ___/\_  __ \/  _ \| __ \_/ __ \ |    |  \_/ __ \ /  ___/  |/ ___\ /    \ 
 |    |     |  | \(  <_> ) \_\ \  ___/ |    `   \  ___/ \___ \|  / /_/  >   |  \\
 |____|     |__|   \____/|___  /\___  >_______  /\___  >____  >__\___  /|___|  /
                             \/     \/        \/     \/     \/  /_____/      \/                             """))

    subparsers = parser.add_subparsers(title='subcommands',
                                       description='valid subcommands',
                                       help='additional help, version: {}'.format(__version__),
                                       dest='action')

    transcript_parse = subparsers.add_parser('transcripts', parents=[probedesign], help='For transcripts ID')
    transcript_parse.add_argument('-faC', '--fastaC', action='store',
                                  type=str,
                                  help='cDNA fasta file, must be modified')
    transcript_parse.set_defaults(func=transcript)

    junction_parse = subparsers.add_parser('junction', parents=[probedesign], help='For splicing junction')
    junction_parse.add_argument('-faG', '--fastaG', action='store', required=True,
                                type=str,
                                help='whole genome fasta file')
    junction_parse.add_argument('-faC', '--fastaC', action='store',
                                type=str,
                                help='cDNA fasta file, must be modified')

    junction_parse.set_defaults(func=junction)

    circ_parse = subparsers.add_parser('circ', parents=[probedesign], help='For circRNA junction')
    circ_parse.add_argument('-faG', '--fastaG', action='store', required=True,
                            type=str,
                            help='whole genome fasta file')
    circ_parse.set_defaults(func=circ)

    # index_parse = subparsers.add_parser('index', parents=[probedesign], help='For index, still in test')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()


def main():
    args = arg()
    args.func(args)


if __name__ == '__main__':
    main()
