import argparse
import sys
import os
import random
import copy
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import readFile,writeFile

parser = argparse.ArgumentParser(description = " generate the simulation ")
parser.add_argument("-i","--input",required=True,help="the input bedpe file.")
parser.add_argument("-o","--outputprefix",required=True,help="the output prefix")
parser.add_argument("-n","--times",required=False,default=10,type=int,help="the simulation times. [default: 10]")
parser.add_argument("-s","--seed",required=False,type=int,help="the random seed")
args = parser.parse_args()

if __name__ == "__main__":
    if args.seed == None:
        seed = random.randint(0,999999999)
        print("you don't give the random seed. Now the random seed is: {0}".format(seed))
    else:
        seed = args.seed
        print("The random seed is: {0}".format(seed))
    random.seed(seed)

    mylist = []
    fin = readFile(args.input)
    ntotal = 0
    for line in fin:
        if line.strip() == "":
            continue
        tmp = line.strip().split()
        mylist.append("\t".join(tmp[0:5]))
        mylist.append("\t".join(tmp[5:10]))
        ntotal += 1
    fin.close()
    
    for ii in range(args.times):
        print("simulation {0} start".format(ii+1))
        fout = writeFile("{0}.sim.{1}.tmp.bedpe".format(args.outputprefix,ii+1))
        newlist = copy.deepcopy(mylist)
        random.shuffle(newlist)
        for index in range(ntotal):
            fout.write("{0}\t{1}\ttmpid{2}\n".format(newlist[2*index],newlist[2*index+1],index))
        fout.close()