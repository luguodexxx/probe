import os
import sys
from collections import defaultdict


def checkdir(dirname):
    """

    :param dirname:
    :return:
    """
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def generateinfo(fasta, targetfile, outputprefix):
    """
    Generate the sequences
    :param fasta:
    :param targetfile:
    :param outputprefix:
    :return:
    """

    targetpool = defaultdict(list)
    fastadict = defaultdict(lambda: defaultdict(list))
    fastalist = []
    status = False
    with open('{}.layerinfo.txt'.format(outputprefix), 'w') as OUT, \
            open(targetfile, 'r') as T, open(fasta, 'r') as FA:

        for gene in T.readlines():
            gene = gene.strip().split('\t')
            targetpool[gene[0]].extend(gene[1:] + [gene[0]])

        for i in FA.readlines():
            if i.startswith('>'):
                i = i.strip().split(' ')
                gid = i[3].split(':')[1].split('.')[0]
                ttype = i[5].split(':')[1]
                tid = i[0].split('|')[0][1:]
                if ttype == 'protein_coding' and gid in targetpool:
                    OUT.write('\t'.join([tid, '\t'.join(targetpool[gid])]) + '\n')
                    fastadict[tid]["header"] = i
                    fastadict[tid]["gid"] = gid
                    status = True
                else:
                    status = False

            else:
                if status:
                    fastadict[tid]["seq"].append(i)
                else:
                    continue
        for k, v in fastadict.items():
            faname = ".".join([k, "fasta"])
            checkdir(v["gid"])
            outprefix = os.path.join(v["gid"], faname)
            fastalist.append(outprefix)
            with open(outprefix, 'w') as FA:
                FA.write(''.join(v["header"]) + '\n')
                FA.write(''.join(v["seq"]))
    return fastalist


def main():
    fasta, targetfile, outputprefix = sys.argv[1:]
    fastalist = generateinfo(fasta, targetfile, outputprefix)


if __name__ == '__main__':
    main()
