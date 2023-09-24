import argparse
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import readFile,writeFile

parser = argparse.ArgumentParser(description = "remove duplicated reads for piared-end reads")
parser.add_argument("-r1","--read1",required=True,help="the input R1 reads")
parser.add_argument("-r2","--read2",required=True,help="the input R2 reads")
parser.add_argument("-o1","--out1",required=True,help="the output R1 read")
parser.add_argument("-o2","--out2",required=True,help="the output R2 read")
args = parser.parse_args()

if __name__ == "__main__":
    R1 = args.read1
    R2 = args.read2
    R1out = args.out1
    R2out = args.out2
    fin1 = readFile(R1)
    fin2 = readFile(R2)
    fout1 = writeFile(R1out)
    fout2 = writeFile(R2out)

    seq = set()
    total = 0
    uniq = 0
    for i,j in zip(fin1,fin2):
        total += 1
        id1 = i
        id2 = j
        seq1 = fin1.readline()
        seq2 = fin2.readline()
        symbol1 = fin1.readline()
        symbol2 = fin2.readline()
        qual1 = fin1.readline()
        qual2 = fin2.readline()

        combineSeq = "{0}#{1}".format(seq1,seq2)
        if combineSeq in seq:
            continue
        else:
            fout1.write("{0}{1}{2}{3}".format(id1,seq1,symbol1,qual1))
            fout2.write("{0}{1}{2}{3}".format(id2,seq2,symbol2,qual2))
            seq.add(combineSeq)
            uniq += 1

    fin1.close()
    fin2.close()
    fout1.close()
    fout2.close()
    sys.stdout.write("The input total reads: {0} and the uniq reads: {1}\n".format(total,uniq))
    sys.stdout.close()