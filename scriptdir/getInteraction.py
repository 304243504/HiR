import argparse
import sys
import os
import random
import copy
import time
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import readFile,writeFile

parser = argparse.ArgumentParser(description = " generate the simulation ")
parser.add_argument("-i","--input",required=True,help="the input bedpe file.")
parser.add_argument("-o","--output",required=True,help="the output file")
parser.add_argument("-l","--genelist",required=True,help="the element file")
args = parser.parse_args()

if __name__ == "__main__":
    coverage = {}
    interact = {}
    fin = readFile(args.input)
    ntotal = 0
    ntotalset = set()
    start = time.time()
    rnadict = {}
    for index,line in enumerate(fin):
        if line.strip() == "":
            continue
        if (index+1) % 1000000 == 0:
            sys.stderr.write('it has load {0} M fragment, use time: {1}s\n'.format((index+1)//1000000,time.time()-start))
            sys.stderr.flush()
            start = time.time()
        tmp = line.strip().split()
        ntotalset.add(tmp[10])
        ntotal += 1
        for ii in tmp[4].split(";"):
            for jj in tmp[9].split(";"):
                if ii != "-":
                    if ii not in coverage:
                        coverage[ii] = []
                    coverage[ii].append(tmp[10])
                    rnadict[ii] = None
                if jj != "-":
                    if jj not in coverage:
                        coverage[jj] = []
                    coverage[jj].append(tmp[10])
                    rnadict[jj] = None
                if ii != "-" and jj != "-":
                    flag = [ii,jj]
                    flag.sort()
                    flagid = "####".join(flag)
                    if flagid not in interact:
                        interact[flagid] = []
                    interact[flagid].append(tmp[10])
    fin.close()
    
    geneinfo = {}
    fin = readFile(args.genelist)
    for line in fin:
        tmp = line.strip().split()
        if tmp[4] == "gene":
            if tmp[3] not in geneinfo:
                geneinfo[tmp[3]] = "{0}\t{1}\t{2}\t{3}\t{4}".format(tmp[0],tmp[1],tmp[2],tmp[3],tmp[5])
            else:
                sys.exit("error! has two same gene, gene id: {0}".format(tmp[3]))
    fin.close()

    sys.stderr.write('print result\n')
    sys.stderr.flush()
    fout = writeFile(args.output)

    for k in rnadict.keys():
        rnadict[k] = "{0}\t{1}".format(len(coverage[k]),len(set(coverage[k])))
    readtotal = len(ntotalset)
    for flagid in interact.keys():
        rna1,rna2 = flagid.split("####")
        fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(geneinfo[rna1],geneinfo[rna2],len(interact[flagid]),len(set(interact[flagid])),rnadict[rna1],rnadict[rna2],ntotal,readtotal))
    fout.close()