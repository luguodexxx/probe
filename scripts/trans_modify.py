#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/24 5:04 PM
__author__ = 'Zhou Ran'

import os
import sys
import re

def main(file):
    prefix = os.path.splitext(file)[0]
    with open(file) as F, open(prefix + '.mod.fa', 'w') as OUT:
        for line in F.readlines():
            if line.startswith(">"):
                gene_id = re.findall('gene:(\w*)', line)
                symbol_id = re.findall('gene_symbol:(\w*)', line)

                line = line.strip().split(" ")
                internal = ":".join(line[2].split(":")[2:])

                OUT.write(' '.join(['|'.join([line[0], symbol_id[0] if len(symbol_id)!=0 else gene_id[0], internal])] + line[1:]) + "\n")
            else:
                OUT.write(line)


if __name__ == '__main__':
    main(sys.argv[1])
