import argparse
import sys
import os
import math
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import readFile,writeFile

parser = argparse.ArgumentParser(description = " generate the simulation ")
parser.add_argument("-i","--input",required=True,help="the input bedpe file.")
parser.add_argument("-ip","--inputprefix",required=True,help="the simulation input prefix")
parser.add_argument("-is","--inputsuffixes",required=True,help="the simulation input suffixes")
parser.add_argument("-n","--times",required=True,type=int,help="the simulation times")
parser.add_argument("-o","--output",required=True,help="the output file")
args = parser.parse_args()

if __name__ == "__main__":
    myhash = {}
    fin = readFile(args.input)
    for line in fin:
        tmp = line.strip().split()
        flagid = "{0}####{1}".format(tmp[3],tmp[8])
        if flagid not in myhash:
            myhash[flagid] = []
        else:
            sys.exit("error! interacion has two same: {0}".format(flagid))
    fin.close()

    for ii in range(args.times):
        fin = readFile("{0}.{1}.{2}".format(args.inputprefix,ii+1,args.inputsuffixes))
        for line in fin:
            tmp = line.strip().split()
            flagid = "{0}####{1}".format(tmp[3],tmp[8])
            if flagid in myhash:
                while len(myhash[flagid]) < ii:
                    myhash[flagid].append(0)
                myhash[flagid].append(int(tmp[10]))
        fin.close()
    
    fin = readFile(args.input)
    fout = writeFile(args.output)
    for line in fin:
        tmp = line.strip().split()
        flagid = "{0}####{1}".format(tmp[3],tmp[8])
        while len(myhash[flagid]) < args.times:
            myhash[flagid].append(0)
        obs = int(tmp[10])
        ngt = 0
        for iiii in myhash[flagid]:
            if iiii >= obs:
                ngt += 1
        VCsqrt = int(tmp[11])/math.sqrt(int(tmp[13])*int(tmp[15]))
        fout.write("{0}\t{1}\t{2}\t{3:.4f}\t{4:.4f}\n".format(line.strip(),ngt,args.times,ngt/args.times,VCsqrt))
    fin.close()
    