#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/24 5:04 PM
__author__ = 'Zhou Ran'

import os
import sys


def main(file):
    prefix = os.path.splitext(file)[0]
    with open(file) as F, open(prefix + '.mod.fa', 'w') as OUT:
        for line in F.readlines():
            if line.startswith(">"):
                line = line.strip().split(" ")
                internal = ":".join(line[2].split(":")[2:])
                genesymbol = line[6].split(":")[1]
                OUT.write(' '.join(['_'.join([line[0], genesymbol, internal])] + line[1:]) + "\n")
            else:
                OUT.write(line)


if __name__ == '__main__':
    main(sys.argv[1])
