import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import readFile,writeFile,regionTree,readSam,getStrand,regionFind,writeSam
import argparse
from collections import OrderedDict

class gene:
    def __init__(self,gid,chrom,start,end,strand):
        self.id = gid
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.rnaTree = {}
    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.chrom,self.start,self.end,self.id,".",self.strand)
    def __repr__(self) -> str:
        pass

class trans:
    def __init__(self,rid,chrom,start,end,strand):
        self.id = rid
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.elementTree = {}
    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.chrom,self.start,self.end,self.id,".",self.strand)
    def __repr__(self) -> str:
        pass

def deal(mylist,geneTree,spliceTree):
    tmpgene = None
    tmprna = None
    for ii in mylist:
        if ii[4] == "gene": # gene
            tmpgene = gene(ii[3],ii[0],ii[1],ii[2],ii[5])
            regionTree(geneTree,ii[0],ii[1],ii[2],tmpgene)
        elif ii[4] == "transcript": # rna
            tmprna = trans(ii[6],ii[0],ii[1],ii[2],ii[5])
            regionTree(tmpgene.rnaTree,ii[0],ii[1],ii[2],tmprna)
        else:
            regionTree(tmprna.elementTree,ii[0],ii[1],ii[2],ii[3:])
            if ii[4] == "intron":
                splicearr = [ii[2]-1,ii[2]]
                splicearr.extend(ii[3:])
                regionTree(spliceTree,ii[0],ii[1],ii[1]+1,splicearr)

def findGene(it,read,readStrand):
    hit = regionFind(it,read.reference_name,read.reference_start,read.reference_end)
    genelist = []
    for ii in hit:
        tmp = ii.value
        if readStrand == tmp.strand:
            genelist.append(ii)
    return genelist

def intersectionRegion(mylist1,mylist2):
    common = []
    for tmpfac in mylist1:
        if tmpfac in mylist2:
            common.append(tmpfac)
    return common

def printTree(mylist):
    for ii in mylist:
        print("{0}\t{1}\t{2}".format(ii.start,ii.end,ii.value))

def checkPairtagRead(read1,read2):
    tmp1 = read1.query_name
    tmp2 = read2.query_name
    if tmp1 != tmp2:
        return True
    return False


parser = argparse.ArgumentParser(description = "split interaction file into intra-gene, inter-gene, without-gene interaction")
parser.add_argument("-i","--input",required=True,help="the input interaction files, comma separated (file1,file2,file3,...  )")
parser.add_argument("-o","--outputprefix",required=True,help="the output prefix")
parser.add_argument("-g","--generegion",required=True,help="the input gene region file")
parser.add_argument("-d","--distance",required=False,type=int,default=500,help="the min distance [default: 500]")
parser.add_argument("-e","--extends",required=False,type=int,default=5,help="the length of exon splice site extend in both of upstream and downstream [default: 5]")
parser.add_argument("-info","--detailedinformation",required=False,choices=["0","1"],default="0",help="print (1) or not print (0) detailed information [default: 0]")
args = parser.parse_args()

if __name__ == "__main__":
    geneTree = {}
    spliceTree = {}
    flag = None
    out = []
    fin = readFile(args.generegion)
    for line in fin:
        tmp = line.strip().split("\t")
        tmp[1] = int(tmp[1])
        tmp[2] = int(tmp[2])
        if len(tmp) > 7:
            tmp[7] = int(tmp[7])
        if flag != None and flag != tmp[3]:
            deal(out,geneTree,spliceTree)
            out = []
        flag = tmp[3]
        out.append(tmp)
    deal(out,geneTree,spliceTree)
    fin.close()
    
    files = args.input.split(",")
    fin = readSam(files[0])
    fsingle = writeSam("{0}.interaction.intraMolecular.Singleton.bam".format(args.outputprefix),fin.header)
    fchim = writeSam("{0}.interaction.intraMolecular.Chimeric.bam".format(args.outputprefix),fin.header)
    flog = writeFile("{0}.interaction.intraMolecular.static.log".format(args.outputprefix))
    if args.detailedinformation == "1":
        finfo = writeFile("{0}.interaction.intraMolecular.info.list".format(args.outputprefix))
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
            readflag = [read1.query_name]
            if strand1 == strand2:
                overlapgene = None
                overlaprna = None
                hit1 = findGene(geneTree,read1,strand1)
                hit2 = findGene(geneTree,read2,strand2)
                overlapgene = intersectionRegion(hit1,hit2)
                if len(overlapgene) == 0:
                    sys.exit("inter and intra class are error! \n{0}\n{1}\n".format(read1.tostring(),read2.tostring()))
                if read1.has_tag("HG") and read2.has_tag("HG"):
                    if len(intersectionRegion(read1.get_tag("HG").split(";"), read2.get_tag("HG").split(";"))) == 0:
                        sys.exit("read's HG tag is not match! \n{0}\n{1}\n".format(read1.tostring(),read2.tostring()))
                    taginfo = read1.get_tag("HG").split(";")
                    if len(taginfo) != len(hit1):
                        printTree(hit1)
                        sys.exit("read's HG tag is not match with IntervalTree1! \n{0}\n".format(read1.tostring()))
                    taginfo = read2.get_tag("HG").split(";")
                    if len(taginfo) != len(hit2):
                        printTree(hit2)
                        sys.exit("read's HG tag is not match with IntervalTree1! \n{0}\n".format(read2.tostring()))
                    for tmpgene in overlapgene:
                        if tmpgene.value.id not in taginfo:
                            printTree(overlapgene)
                            sys.exit("read's HG tag is not match with IntervalTree2! \n{0}\n".format(read2.tostring()))
                for tmpgene in overlapgene:
                    rnahit1 = findGene(tmpgene.value.rnaTree,read1,strand1)
                    rnahit2 = findGene(tmpgene.value.rnaTree,read2,strand2)
                    overlaprna = intersectionRegion(rnahit1,rnahit2)
                if strand1 == "+":
                    if read1.reference_start - read2.reference_end < args.distance:
                        if read1.reference_start - read2.reference_end > -100: # singleton
                            read1.set_tag("PT","inpreRNAshorterthan{0}".format(args.distance))  # pairTag type
                            read2.set_tag("PT","inpreRNAshorterthan{0}".format(args.distance))
                            pairType = "inpreRNAshorterthan{0}".format(args.distance)
                            readflag.append(pairType)
                            readflag.append("{0}".format(read1.reference_start - read2.reference_end))
                            fsingle.write(read1)
                            fsingle.write(read2)
                        else:
                            read1.set_tag("PT","locationchange")
                            read2.set_tag("PT","locationchange")
                            pairType = "locationchange"
                            readflag.append(pairType)
                            readflag.append("{0}".format(read1.reference_start - read2.reference_end))
                            fchim.write(read1)
                            fchim.write(read2)
                    elif len(overlaprna) == 0:
                        read1.set_tag("PT","NoRNA")
                        read2.set_tag("PT","NoRNA")
                        pairType = "NoRNA"
                        readflag.append("{0}".format(read1.reference_start - read2.reference_end))
                        readflag.append(pairType)
                        fchim.write(read1)
                        fchim.write(read2)
                    else:
                        fragmentlen = []
                        for tmprna in overlaprna:
                            regionhit1 = regionFind(tmprna.value.elementTree,read1.reference_name,read1.reference_start,read1.reference_start+1)
                            regionhit2 = regionFind(tmprna.value.elementTree,read2.reference_name,read2.reference_end-1,read2.reference_end)
                            if len(regionhit1) != 1 or len(regionhit2) != 1:
                                sys.exit("reads has multiple region! \n{0}\n{1}\n{2}\n{3}\n".format(read1.tostring(),read2.tostring(),regionhit1,regionhit2))
                            if regionhit1[0].value[1] == "intron" or regionhit2[0].value[1] == "intron":
                                RNAdis = read1.reference_start - read2.reference_end
                                fragmentlen.append([RNAdis,tmprna.value.id,"intron"])
                            else:
                                tss1 = regionhit1[0].value[4] + read1.reference_start - regionhit1[0].start
                                tss2 = regionhit2[0].value[4] + read2.reference_end - regionhit2[0].start
                                RNAdis = tss1 - tss2
                                fragmentlen.append([RNAdis,tmprna.value.id,"exon"])
                        fragmentlen.sort(key=lambda x:x[0])
                        if fragmentlen[0][0] < args.distance: # singleton
                            read1.set_tag("PT","inMatureRNAshorterthan{0}in{1}".format(args.distance,fragmentlen[0][2]))
                            read2.set_tag("PT","inMatureRNAshorterthan{0}in{1}".format(args.distance,fragmentlen[0][2]))
                            pairType = "inMatureRNAshorterthan{0}in{1}".format(args.distance,fragmentlen[0][2])
                            readflag.append(pairType)
                            readflag.append("{0}".format(read1.reference_start - read2.reference_end))
                            for iiii in fragmentlen:
                                readflag.append("{0}#{1}#{2}".format(iiii[0],iiii[1],iiii[2]))
                            fsingle.write(read1)
                            fsingle.write(read2)
                        else:
                            read1.set_tag("PT","inMatureRNAlongerthan{0}in{1}".format(args.distance,fragmentlen[0][2]))
                            read2.set_tag("PT","inMatureRNAlongerthan{0}in{1}".format(args.distance,fragmentlen[0][2]))
                            pairType = "inMatureRNAlongerthan{0}in{1}".format(args.distance,fragmentlen[0][2])
                            readflag.append(pairType)
                            readflag.append("{0}".format(read1.reference_start - read2.reference_end))
                            for iiii in fragmentlen:
                                readflag.append("{0}#{1}#{2}".format(iiii[0],iiii[1],iiii[2]))
                            fchim.write(read1)
                            fchim.write(read2)
                else:
                    if read2.reference_start - read1.reference_end < args.distance:
                        if read2.reference_start - read1.reference_end > -100: # singleton
                            read1.set_tag("PT","inpreRNAshorterthan{0}".format(args.distance))  # pairTag type
                            read2.set_tag("PT","inpreRNAshorterthan{0}".format(args.distance))
                            pairType = "inpreRNAshorterthan{0}".format(args.distance)
                            readflag.append(pairType)
                            readflag.append("{0}".format(read2.reference_start - read1.reference_end))
                            fsingle.write(read1)
                            fsingle.write(read2)
                        else:
                            read1.set_tag("PT","locationchange")
                            read2.set_tag("PT","locationchange")
                            pairType = "locationchange"
                            readflag.append(pairType)
                            readflag.append("{0}".format(read2.reference_start - read1.reference_end))
                            fchim.write(read1)
                            fchim.write(read2)
                    elif len(overlaprna) == 0:
                        read1.set_tag("PT","NoRNA")
                        read2.set_tag("PT","NoRNA")
                        pairType = "NoRNA"
                        readflag.append("{0}".format(read2.reference_start - read1.reference_end))
                        readflag.append(pairType)
                        fchim.write(read1)
                        fchim.write(read2)
                    else:
                        fragmentlen = []
                        for tmprna in overlaprna:
                            regionhit1 = regionFind(tmprna.value.elementTree,read1.reference_name,read1.reference_end-1,read1.reference_end)
                            regionhit2 = regionFind(tmprna.value.elementTree,read2.reference_name,read2.reference_start,read2.reference_start+1)
                            if len(regionhit1) != 1 or len(regionhit2) != 1:
                                sys.exit("reads has multiple region! \n{0}\n{1}\n{2}\n{3}\n".format(read1.tostring(),read2.tostring(),regionhit1,regionhit2))
                            if regionhit1[0].value[1] == "intron" or regionhit2[0].value[1] == "intron":
                                RNAdis = read2.reference_start - read1.reference_end
                                fragmentlen.append([RNAdis,tmprna.value.id,"intron"])
                            else:
                                tss1 = regionhit1[0].value[4] + regionhit1[0].end - read1.reference_end
                                tss2 = regionhit2[0].value[4] + regionhit2[0].end - read2.reference_start
                                RNAdis = tss1 - tss2
                                fragmentlen.append([RNAdis,tmprna.value.id,"exon"])
                        fragmentlen.sort(key=lambda x:x[0])
                        if fragmentlen[0][0] < args.distance: # singleton
                            read1.set_tag("PT","inMatureRNAshorterthan{0}in{1}".format(args.distance,fragmentlen[0][2]))
                            read2.set_tag("PT","inMatureRNAshorterthan{0}in{1}".format(args.distance,fragmentlen[0][2]))
                            pairType = "inMatureRNAshorterthan{0}in{1}".format(args.distance,fragmentlen[0][2])
                            readflag.append(pairType)
                            readflag.append("{0}".format(read2.reference_start - read1.reference_end))
                            for iiii in fragmentlen:
                                readflag.append("{0}#{1}#{2}".format(iiii[0],iiii[1],iiii[2]))
                            fsingle.write(read1)
                            fsingle.write(read2)
                        else:
                            read1.set_tag("PT","inMatureRNAlongerthan{0}in{1}".format(args.distance,fragmentlen[0][2]))
                            read2.set_tag("PT","inMatureRNAlongerthan{0}in{1}".format(args.distance,fragmentlen[0][2]))
                            pairType = "inMatureRNAlongerthan{0}in{1}".format(args.distance,fragmentlen[0][2])
                            readflag.append(pairType)
                            readflag.append("{0}".format(read2.reference_start - read1.reference_end))
                            for iiii in fragmentlen:
                                readflag.append("{0}#{1}#{2}".format(iiii[0],iiii[1],iiii[2]))
                            fchim.write(read1)
                            fchim.write(read2)
            else:
                read1.set_tag("PT","diffstrand")
                read2.set_tag("PT","diffstrand")
                pairType = "diffstrand"
                readflag.append(pairType)
                fchim.write(read1)
                fchim.write(read2)
            if args.detailedinformation == "1":
                finfo.write("{0}\n".format("\t".join(readflag)))
            if pairType not in myhash:
                myhash[pairType] = 0
            myhash[pairType] += 1
        fin.close()
    fsingle.close()
    fchim.close()
    if args.detailedinformation == "1":
        finfo.close()
    for ii in sorted(myhash.keys()):
        flog.write("{0}\t{1}\n".format(ii,myhash[ii]))
    flog.close()