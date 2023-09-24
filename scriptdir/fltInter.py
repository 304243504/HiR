import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import writeFile,readSam,getStrand,writeSam,readFile,regionTree
import argparse
from bx.intervals.intersection import Intersecter, Interval

def buildTree(resFrag,chrom,start,end,strand,flag):
    if chrom not in resFrag:
        resFrag[chrom] = Intersecter()
    tmpi = Interval(start,end,flag)
    tmpi.strand = strand
    resFrag[chrom].add_interval(tmpi)

def regionFind(it,chrom,start,end,strand):
    if chrom not in it:
        return []
    tmpr = it[chrom].find(start,end)
    resultarr = []
    for ii in tmpr:
        if ii.strand == strand:
            resultarr.append(ii)
    return resultarr

def regionFind2(it,chrom,start,end):
    if chrom not in it:
        return []
    return it[chrom].find(start,end)

def overlapflag(itl,chrom,start,end,strand):
    if itl.value[0] == chrom and itl.value[3] == strand:
        return itl.value[2] > start and itl.value[1] < end
    return False

def checkPairtagRead(read1,read2):
    tmp1 = read1.query_name
    tmp2 = read2.query_name
    if tmp1 != tmp2:
        return True
    return False

parser = argparse.ArgumentParser(description = "split interaction file into intra-gene, inter-gene, without-gene interaction")
parser.add_argument("-i","--input",required=True,help="the input interaction files, comma separated (file1,file2,file3,...  )")
parser.add_argument("-r","--rnaseq",required=True,help="the input RNA-seq or totalRNA-seq interaction files as control, comma separated (file1,file2,file3,...  )")
parser.add_argument("-c","--cutoff",required=False,type=int,default=4,help=" for RNA-Seq the min chimeric number [default: 4]")
parser.add_argument("-g","--generegion",required=True,help="the input gene region file")
parser.add_argument("-d","--distance",required=False,type=int,default=500,help="the min distance [default: 500]")
parser.add_argument("-o","--outputprefix",required=True,help="the output prefix")
args = parser.parse_args()

if __name__ == "__main__":
    spliceTree = {}
    fin = readFile(args.generegion)
    for line in fin:
        tmp = line.strip().split("\t")
        tmp[1] = int(tmp[1])
        tmp[2] = int(tmp[2])
        if tmp[4] == "intron":
            splicearr = [tmp[2]-1,tmp[2]] # splice site
            splicearr.extend(tmp[3:])
            regionTree(spliceTree,tmp[0],tmp[1],tmp[1]+1,splicearr) # splice site
    fin.close()

    it = {} # interval tree
    for tmpfile in args.rnaseq.split(","):
        fin = readFile(tmpfile)
        for line in fin:
            tmp = line.strip().split()
            if int(tmp[8]) < args.cutoff:
                continue
            record = [[tmp[0],int(tmp[1]),int(tmp[2]),tmp[6]],[tmp[3],int(tmp[4]),int(tmp[5]),tmp[7]]]
            record[1].append(tmp[-1].split(";")[0])
            buildTree(it,record[0][0],record[0][1],record[0][2],record[0][3],record[1])
        fin.close()

    files = args.input.split(",")
    fin = readSam(files[0])
    fdrop = writeSam("{0}.interaction.interMolecular.Drop.bam".format(args.outputprefix),fin.header)
    fchim = writeSam("{0}.interaction.interMolecular.Chimeric.bam".format(args.outputprefix),fin.header)
    fsingle = writeSam("{0}.interaction.interMolecular.Singleton.bam".format(args.outputprefix),fin.header)
    flog = writeFile("{0}.interaction.interMolecular.static.log".format(args.outputprefix))
    fin.close()
    myhash = {}

    for tmpfile in files:
        fin = readSam(tmpfile)
        for read1 in fin:
            read2 = next(fin)
            if checkPairtagRead(read1,read2):
                sys.stderr.write("two reads are not match \n{0}\n{1}\n".format(read1.tostring(),read2.tostring()))
            if read1.is_reverse:
                strand1 = "-"
            else:
                strand1 = "+"
            if read2.is_reverse:
                strand2 = "-"
            else:
                strand2 = "+"
            pairType = None
            is_remove = True
            if is_remove:
                if read1.is_reverse:
                    ss = read1.reference_start
                else:
                    ss = read1.reference_end - 1
                        
                if read2.is_reverse:
                    ee = read2.reference_end - 1
                else:
                    ee = read2.reference_start
              record = [[read1.reference_name,ss-5,ss+6,strand1],[read2.reference_name,ee-5,ee+6,strand2]] # extends 5bp
                record.sort(key=lambda x:(x[0],x[1],x[2],x[3]))
                hitlist = regionFind(it,record[0][0],record[0][1],record[0][2],record[0][3])

                if len(hitlist) == 0:
                    fchim.write(read1)
                    fchim.write(read2)
                    pairType = "Chimeric"
                else:
                    flagt = False
                    for iii in hitlist:
                        if overlapflag(iii,record[1][0],record[1][1],record[1][2],record[1][3]):
                            flagt = True
                            read1.set_tag("IO",iii.value[4]) # interMolecular overlap
                            read2.set_tag("IO",iii.value[4])
                            break
                    if flagt:
                        fdrop.write(read1)
                        fdrop.write(read2)
                        pairType = "Drop"
                    else:
                        fchim.write(read1)
                        fchim.write(read2)
                        pairType = "Chimeric"
            if pairType not in myhash:
                myhash[pairType] = 0
            myhash[pairType] += 1
        fin.close()
    fdrop.close()
    fchim.close()
    fsingle.close()
    for ii in sorted(myhash.keys()):
        flog.write("{0}\t{1}\n".format(ii,myhash[ii]))
    flog.close()
