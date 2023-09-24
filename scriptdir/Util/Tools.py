"""
20220719: add Haplotype
"""
### Operating common files  ###
import gzip
def readFile(infile):
    """
    infile: input file
    return: file handle
    """
    if infile.endswith((".gz","gzip",".GZ",".GZIP")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin

def writeFile(outfile):
    """
    outfile: output file
    return: file handle
    """
    if outfile.endswith((".gz","gzip",".GZ",".GZIP")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

### Operating sam/bam files  ###
import pysam
def readSam(insamfile):
    """
    insamfile: input sam/bam file
    return: file handle
    """
    if insamfile.endswith((".bam",".sam.gz")):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam"):
        insam = pysam.AlignmentFile(insamfile,'r')
    else:
        raise ValueError("the input sam/bam file is not end with sam or bam!")
    return insam

def writeSam(outsamfile,header):
    """
    outsamfile: output sam/bam file
    header: the sam/bam file's header(chromosome information, created by insam.handle)
    return: file handle
    """
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

def addGroup(mylist,tag,taginfo):
    for read in mylist:
        read.set_tag(tag,taginfo)

numberset = set(["0","1","2","3","4","5","6","7","8","9"])
def MDstr2MDarr(mystr):
    MDarrtmp = []
    flag = ""
    for i in mystr:
        if i not in numberset:
            if flag != "":
                MDarrtmp.append(int(flag))
            MDarrtmp.append(i)
            flag = ""
        else:
            flag = flag+i
    if flag != "":
        MDarrtmp.append(int(flag))
    return MDarrtmp

def MDarr2MDstr(mylist):
    tmpstr = ""
    for i in mylist:
        tmpstr = "{0}{1}".format(tmpstr,i)
    return tmpstr

import sys
def getReadPos(read):
    arr = read.cigartuples
    if read.is_reverse:
        arr = arr[::-1]
    start = 0
    end = 0
    for index,(i,j) in enumerate(arr):
        if index == 0:
            if i==4 or i==5: # S H soft/hard clip
                start += j
                end += j
            elif i == 1: # I insertion
                # grep E00513:271:HGKL2CCX2:3:2105:13210:17166 align/RIC-013/RIC-013.genome.read2.futher_by_bwa.sam
                # E00513:271:HGKL2CCX2:3:2105:13210:17166	0	Chr06	26548988	23	3I84M18S	*	0	0	GTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGAGGTGGCGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGCGGCGGTGGTGGTG	AAFFFJJFJJFJJFJJJJJAJJJJJFFJFJJFJJFJJ7FJ-FJFJJAJJ-77-AJ-<J<JJ-<F7FA7FFJJJAFF7<F-7F-7F-A<-FJ<JJ-7A-77AF--7	NM:i:8	MD:Z:33A18A7C8A6A7	AS:i:50	XS:i:38	SA:Z:Chr05,1524008,+,73S32M,0,1;
                # E00513:271:HGKL2CCX2:3:2105:13210:17166	2048	Chr05	1524008	0	73H32M	*	0	0	TGGTGGTGGTGGTGGTGGCGGCGGTGGTGGTG	AFF7<F-7F-7F-A<-FJ<JJ-7A-77AF--7	NM:i:1	MD:Z:29A2	AS:i:29	XS:i:28	SA:Z:Chr06,26548988,+,3I84M18S,23,8;
                # sys.exit("Impossible1! The cigar: {0}\t{1}".format(read.cigarstring,read.to_string()))
                end += j
                sys.stderr.write("{0} have strange mapping fragment\n".format(read.query_name))
            elif i == 2: # D deletion
                sys.exit("Impossible2! The cigar: {0}\n{1}".format(read.cigarstring,read.to_string()))
            elif i == 0: # M match
                start += 0
                end += j
            else:
                sys.exit("Impossible3! The cigar: {0}\n{1}".format(read.cigarstring,read.to_string()))
        else:
            if i==4 or i==5: # S H soft/hard clip
                pass
            elif i == 1: # I insertion
                end += j
            elif i == 2: # D deletion
                pass
            elif i==3: # N skip region
                pass
            elif i == 0: # M match
                end += j
            else:
                sys.exit("Impossible3! The cigar: {0}\n{1}".format(read.cigarstring,read.to_string()))
    return start,end

def getReadPosPro(read):
    arr = read.cigartuples
    if read.is_reverse:
        arr = arr[::-1]
    start = 0
    end = 0
    newcigar = []
    for index,(i,j) in enumerate(arr):
        if index == 0:
            if i==4 or i==5: # S H soft/hard clip
                start += j
                end += j
            elif i == 1: # I insertion
                end += j
                newcigar.append((i,j))
                sys.stderr.write("{0} have strange mapping fragment\n".format(read.query_name))
            elif i == 2: # D deletion
                newcigar.append((i,j))
                sys.exit("Impossible2! The cigar: {0}\n{1}".format(read.cigarstring,read.to_string()))
            elif i == 0: # M match
                start += 0
                end += j
                newcigar.append((i,j))
            else:
                sys.exit("Impossible3! The cigar: {0}\n{1}".format(read.cigarstring,read.to_string()))
        else:
            if i==4 or i==5: # S H soft/hard clip
                pass
            elif i == 1: # I insertion
                end += j
                newcigar.append((i,j))
            elif i == 2: # D deletion
                newcigar.append((i,j))
            elif i==3: # N skip region
                newcigar.append((i,j))
            elif i == 0: # M match
                end += j
                newcigar.append((i,j))
            else:
                sys.exit("Impossible3! The cigar: {0}\n{1}".format(read.cigarstring,read.to_string()))
    if read.is_reverse:
        newcigar = newcigar[::-1]
    return start,end,newcigar

def coutFragmentDis(s1,e1,s2,e2): # >0: distance;  <0: overlap
    min = s1
    if min < s2:
        min = s2
    max = e1
    if max > e2:
        max = e2
    return min-max

def SEsam2fastq(read):
    return "@{0}\n{1}\n+\n{2}\n".format(read.query_name,read.get_forward_sequence(),pysam.qualities_to_qualitystring(read.get_forward_qualities()))

def qual2str(read):
    return pysam.qualities_to_qualitystring(read.get_forward_qualities())

from collections import OrderedDict
def parseGFFString(mystr):
    """
    mystr: gtf attribution information string
    return: a dictionary
    """
    d = OrderedDict()
    if mystr.endswith(";"):
        mystr = mystr[:-1]
    for i in mystr.strip().split(";"):
        tmp = i.strip().split("=")
        if tmp == "":
            continue
        d[tmp[0]] = tmp[1]
    return d

def rc(sequence):
    seq = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = seq.translate(trantab)
    return string

def checkChimeric(mylist):
    begin1,end1 = getReadPos(mylist[0])
    begin2,end2 = getReadPos(mylist[1])
    dis = coutFragmentDis(begin1,end1,begin2,end2)
    maplen = end1-begin1+end2-begin2
    if dis < 0:
        maplen += dis
    beforeS = begin1
    afterS = mylist[0].query_length - end2
    return mylist[0].query_length,beforeS,afterS,maplen,dis

def uniqstar(mylist):
    if len(mylist) > 1:
        return "multimap"
    elif len(mylist) == 1:
        if mylist[0].is_unmapped:
            return "unmap"
        else:
            return "uniqmap"
    else:
        for ii in mylist:
            print(ii)
        sys.exit("error")

def uniqstarsingleread(read):
    if read.is_unmapped:
        return "unmap"
    else:
        if read.get_tag("NH") == 1:
            return "uniqmap"
        else:
            return "multimap"

def mapstatbwa(read):
    if read.is_unmapped:
        return "unmap"
    elif read.get_tag("AS") != read.get_tag("XS"): # uniq reads
        return "uniqmap"
    else:
        return "multimap"

def toDict(header):
    mydict = {}
    for i in range(header.nreferences):
        mydict[header.get_reference_name(i)] = header.get_reference_length(header.get_reference_name(i))
    return mydict

def starAlignLength(mylist):
    begin1,end1 = getReadPos(mylist[0])
    maplen = end1-begin1
    beforeS = begin1
    afterS = mylist[0].query_length - end1
    return mylist[0].query_length,beforeS,afterS,maplen,0

def checkCSB(chimericreads,starreads,bwareads):
    flags = False
    for ii in chimericreads:
        if flags == True:
            break
        for jj in starreads:
            if jj.is_unmapped:
                continue
            if ii.reference_name == jj.reference_name:
                if ii.reference_start - jj.reference_start > -100 and ii.reference_start - jj.reference_start < 100:
                    flags = True
                    break
                if ii.reference_end - jj.reference_end > -100 and ii.reference_end - jj.reference_end < 100:
                    flags = True
                    break
    if starreads[0].is_unmapped:
        aa = ["2","##"]
    else:
        if flags == True:
            aa = ["1","##"]
        else:
            aa = ["0","##"]
    flags = False
    for ii in chimericreads:
        if flags == True:
            break
        for jj in bwareads:
            if jj.is_unmapped:
                continue
            if ii.reference_name == jj.reference_name:
                if ii.reference_start - jj.reference_start > -100 and ii.reference_start - jj.reference_start < 100:
                    flags = True
                    break
                if ii.reference_end - jj.reference_end > -100 and ii.reference_end - jj.reference_end < 100:
                    flags = True
                    break

    if bwareads[0].is_unmapped:
        aa[1] = "2"
    else:
        if flags == True:
            aa[1] = "1"
        else:
            aa[1] = "0"
    return aa

def checkCS(chimericreads,starreads):
    flags = False
    for ii in chimericreads:
        if flags == True:
            break
        for jj in starreads:
            if jj.is_unmapped:
                continue
            if ii.reference_name == jj.reference_name:
                if ii.reference_start - jj.reference_start > -100 and ii.reference_start - jj.reference_start < 100:
                    flags = True
                    break
                if ii.reference_end - jj.reference_end > -100 and ii.reference_end - jj.reference_end < 100:
                    flags = True
                    break
    if starreads[0].is_unmapped:
        aa = ["2","##"]
    else:
        if flags == True:
            aa = ["1","##"]
        else:
            aa = ["0","##"]
    return aa

def checkSB(starreads,bwareads):
    flags1 = False
    flags2 = False
    if starreads[0].is_unmapped and bwareads[0].is_unmapped:
        return ["2","2"]
    if starreads[0].is_unmapped and not bwareads[0].is_unmapped:
        return ["2","0"]
    if not starreads[0].is_unmapped and bwareads[0].is_unmapped:
        return ["0","2"]
    for ii in starreads:
        if flags1 == True and flags2 == True:
            break
        for jj in bwareads:
            if ii.reference_name == jj.reference_name:
                if ii.reference_start - jj.reference_start > -100 and ii.reference_start - jj.reference_start < 100:
                    flags1 = True
                if ii.reference_end - jj.reference_end > -100 and ii.reference_end - jj.reference_end < 100:
                    flags2 = True
    if flags1 == True:
        aa = ["1","##"]
    else:
        aa = ["0","##"]
    if flags2 == True:
        aa[1] = "1"
    else:
        aa[1] = "0"
    return aa

def checkAlignUniq(mylist):
    begin1,end1 = getReadPos(mylist[0])
    maplen = end1-begin1
    beforeS = begin1
    afterS = mylist[0].query_length - end1
    return mylist[0].query_length,beforeS,afterS,maplen,0

def checkStarAlign(reads):
    if len(reads) > 1: # multi mapping
        index = 0
        for index in range(len(reads)):
            if not reads[index].is_secondary:
                break
        reads[0],reads[index] = reads[index], reads[0]
        if reads[0].is_secondary:
            for ii in reads:
                print(ii)
            sys.exit("please use the raw star align file")
        readlen,beforeS,afterS,maplen,dis = checkAlignUniq(reads[:1])
        if beforeS < 15 and afterS < 15: # don't need realign
            # return "multimap#full"
            return ["multimap#full", beforeS, afterS]
        else:
            # return "multimap#notfull"
            return ["multimap#notfull", beforeS, afterS]
    else:
        if reads[0].is_unmapped: # unmap mapping, 
            # return "unmap#full"
            return ["unmap#full", 0, 0]
        else: # uniq mapping
            if reads[0].get_tag("NH") != 1:
                sys.exit("please check STAR's --outFilterMultimapNmax parameters ")
            readlen,beforeS,afterS,maplen,dis = checkAlignUniq(reads)
            if beforeS < 15 and afterS < 15: # don't need realign
                # return "uniqmap#full"
                return ["uniqmap#full", beforeS, afterS]
            else:
                # return "uniqmap#notfull"
                return ["uniqmap#notfull", beforeS, afterS]

def checkBwaAlign(reads):
    if len(reads) > 1: # splice reads
        blocks = bwaRead2Blocks(reads)
        beforeS = blocks[0].start
        afterS = len(blocks[0].raw_seq) - blocks[-1].end
        if beforeS < 15 and afterS < 15: # don't need realign
            # return "uniqormultimap#full"
            return ["uniqormultimap#full", beforeS, afterS]
        else:
            # return "uniqormultimap#notfull"
            return ["uniqormultimap#notfull", beforeS, afterS]
    else:
        if reads[0].is_unmapped: # unmap mapping, 
            # return "unmap#full"
            return ["unmap#full", 0, 0]
        else: # uniq mapping
            readlen,beforeS,afterS,maplen,dis = checkAlignUniq(reads[:1])
            if reads[0].get_tag("AS") != reads[0].get_tag("XS"):
                if beforeS < 15 and afterS < 15: # don't need realign
                    # return "uniqmap#full"
                    return ["uniqmap#full", beforeS, afterS]
                else:
                    # return "uniqmap#notfull"
                    return ["uniqmap#notfull", beforeS, afterS]
            else:
                if beforeS < 15 and afterS < 15: # don't need realign
                    # return "multimap#full"
                    return ["multimap#full", beforeS, afterS]
                else:
                    # return "multimap#notfull"
                    return ["multimap#notfull", beforeS, afterS]

class block:
    def __init__(self,read,qstart,qend,rstart,rend,newcigar,mapstat="uniqmap",raw_seq=None,raw_qual=None,raw_cigar=None,MD=None,raw_MD=None):
        self.read = read
        self.start = qstart
        self.end = qend
        self.rstart = rstart
        self.rend = rend
        self.newcigar = newcigar
        self.mapstat = mapstat
        self.raw_seq = raw_seq
        self.raw_qual = raw_qual
        self.raw_cigar = raw_cigar
        self.MD = MD
        self.raw_MD = raw_MD
    def __str__(self):
        return "{0}\t{1}\t{2}".format(self.start,self.end,self.read.to_string())
    def __repr__(self):
        return self.__str__()

def starRead2Blocks(read):
    blocks = []
    qstart = 0
    qend = 0
    mycigar = []
    arr = read.cigartuples
    raw_seq = read.get_forward_sequence()
    raw_qual = pysam.qualities_to_qualitystring(read.get_forward_qualities())
    raw_cigar = read.cigarstring
    MDarr = MDstr2MDarr(read.get_tag("MD"))
    # print(read)
    if read.is_reverse:
        arr = arr[::-1]
        MDrest = MDarr
        MDcur = ""
        cstart = read.reference_end
        cend = read.reference_end
        for index,(i,j) in enumerate(arr):
            if i < 3:
                mycigar.append((i,j))
                if i == 0: # M
                    qend += j
                    cstart -= j
                    MDpos = 0
                    while MDpos != j:
                        # print(MDrest,MDpos,MDcur,j)
                        if MDrest[-1].__class__ == int:
                            if MDrest[-1] == 0:
                                # sys.stderr.write("{0}\n".format(read.query_name))
                                # sys.stderr.flush()
                                MDcur = "{0}{1}".format(MDrest.pop(),MDcur)
                                continue
                            if MDpos + MDrest[-1] <= j:
                                tmptmptmp = MDrest.pop()
                                MDpos += tmptmptmp
                                MDcurarr = MDstr2MDarr(MDcur)
                                if MDcur == "":
                                    MDcur = "{0}{1}".format(tmptmptmp,MDcur)
                                elif MDcurarr[0].__class__ == int:
                                    MDcurarr[0] += tmptmptmp
                                    MDcur = MDarr2MDstr(MDcurarr)
                                else:
                                    MDcur = "{0}{1}".format(tmptmptmp,MDcur)
                            else:
                                tmptmptmp = j-MDpos
                                MDpos += tmptmptmp
                                MDcurarr = MDstr2MDarr(MDcur)
                                # print(MDcurarr)
                                if MDcur == "":
                                    MDcur = "{0}{1}".format(tmptmptmp,MDcur)
                                elif MDcurarr[0].__class__ == int:
                                    MDcurarr[0] += tmptmptmp
                                    MDcur = MDarr2MDstr(MDcurarr)
                                else:
                                    MDcur = "{0}{1}".format(tmptmptmp,MDcur)
                                MDrest[-1] = MDrest[-1] - tmptmptmp
                                # print(MDcur,MDrest)
                        else:
                            tmptmptmp = MDrest.pop()
                            MDpos += 1
                            MDcur = "{0}{1}".format(tmptmptmp,MDcur)
                elif i == 1: # I insertion
                    qend += j
                    if index == 0:
                        sys.stderr.write("{0} have strange mapping fragment\n".format(read.query_name))
                        sys.stderr.flush()
                elif i == 2: # D deletion
                    cstart -= j
                    if MDrest[-1-j] != "^":
                        sys.stderr.write("{0} has the error MD flag".format(read.query_name))
                        sys.stderr.flush()
                        sys.exit("please check! mark15")
                    MDpos = 0
                    while MDpos != j:
                        if MDrest[-1].__class__ == int:
                            sys.stderr.write("{0} has the error MD flag".format(read.query_name))
                            sys.stderr.flush()
                            sys.exit("please check! mark16")
                        else:
                            tmptmptmp = MDrest.pop()
                            MDpos += 1
                            MDcur = "{0}{1}".format(tmptmptmp,MDcur)
                    MDcur = "{0}{1}".format(MDrest.pop(),MDcur)
                    if index == 0:
                        sys.stderr.write("{0} have strange mapping fragment\n".format(read.query_name))
                        sys.stderr.flush()
            elif i == 3: # N skip region
                blocks.append(block(read,qstart,qend,cstart,cend,mycigar[::-1],uniqstarsingleread(read),raw_seq,raw_qual,raw_cigar,MDcur,read.get_tag("MD")))
                MDcur = ""
                mycigar = []
                qstart = qend
                cend = cstart - j
                cstart = cend
            elif i == 4 or i == 5: # S H soft/hard clip
                if index == 0:
                    qstart += j
                    qend += j
            else:
                sys.exit("Impossible3! The cigar error: {0}\n{1}".format(read.cigarstring,read.to_string()))
        if not (len(MDrest) == 0 or (len(MDrest) == 1 and MDrest[0] == 0)):
            print(MDrest)
            print(MDcur)
            sys.stderr.write("{0} has the error MD flag".format(read.query_name))
            sys.stderr.flush()
            sys.exit("please check! mark17")
        blocks.append(block(read,qstart,qend,cstart,cend,mycigar[::-1],uniqstarsingleread(read),raw_seq,raw_qual,raw_cigar,MDcur,read.get_tag("MD")))
        MDcur = ""
    else:
        MDrest = MDarr[::-1]
        MDcur = ""
        cstart = read.reference_start
        cend = read.reference_start
        for index,(i,j) in enumerate(arr):
            if i < 3:
                mycigar.append((i,j))
                if i == 0: # M
                    qend += j
                    cend += j
                    MDpos = 0
                    while MDpos != j:
                        if MDrest[-1].__class__ == int:
                            if MDrest[-1] == 0:
                                # sys.stderr.write("{0}\n".format(read.query_name))
                                # sys.stderr.flush()
                                MDcur = "{0}{1}".format(MDcur,MDrest.pop())
                                continue
                            if MDpos + MDrest[-1] <= j:
                                tmptmptmp = MDrest.pop()
                                MDpos += tmptmptmp
                                MDcurarr = MDstr2MDarr(MDcur)
                                if MDcur == "":
                                    MDcur = "{0}{1}".format(MDcur,tmptmptmp)
                                elif MDcurarr[-1].__class__ == int:
                                    MDcurarr[-1] += tmptmptmp
                                    MDcur = MDarr2MDstr(MDcurarr)
                                else:
                                    MDcur = "{0}{1}".format(MDcur,tmptmptmp)
                            else:
                                tmptmptmp = j-MDpos
                                MDpos += tmptmptmp
                                MDcurarr = MDstr2MDarr(MDcur)
                                if MDcur == "":
                                    MDcur = "{0}{1}".format(MDcur,tmptmptmp)
                                elif MDcurarr[-1].__class__ == int:
                                    MDcurarr[-1] += tmptmptmp
                                    MDcur = MDarr2MDstr(MDcurarr)
                                else:
                                    MDcur = "{0}{1}".format(MDcur,tmptmptmp)
                                MDrest[-1] = MDrest[-1] - tmptmptmp
                        else:
                            tmptmptmp = MDrest.pop()
                            MDpos += 1
                            MDcur = "{0}{1}".format(MDcur,tmptmptmp)
                elif i == 1: # I insertion
                    qend += j
                    if index == 0:
                        sys.stderr.write("{0} have strange mapping fragment\n".format(read.query_name))
                        sys.stderr.flush()
                elif i == 2: # D deletion
                    cend += j
                    if MDrest[-1] != "^":
                        sys.stderr.write("{0} has the error MD flag".format(read.query_name))
                        sys.stderr.flush()
                        sys.exit("please check! mark18")
                    MDcur = "{0}{1}".format(MDcur,MDrest.pop())
                    MDpos = 0
                    while MDpos != j:
                        if MDrest[-1].__class__ == int:
                            sys.stderr.write("{0} has the error MD flag".format(read.query_name))
                            sys.stderr.flush()
                            sys.exit("please check! mark19")
                        else:
                            tmptmptmp = MDrest.pop()
                            MDpos += 1
                            MDcur = "{0}{1}".format(MDcur,tmptmptmp)
                    if index == 0:
                        sys.stderr.write("{0} have strange mapping fragment\n".format(read.query_name))
                        sys.stderr.flush()
            elif i == 3: # N skip region
                blocks.append(block(read,qstart,qend,cstart,cend,mycigar,uniqstarsingleread(read),raw_seq,raw_qual,raw_cigar,MDcur,read.get_tag("MD")))
                MDcur = ""
                mycigar = []
                qstart = qend
                cstart = cend + j
                cend = cstart
            elif i == 4 or i == 5: # S H soft/hard clip
                if index == 0:
                    qstart += j
                    qend += j
            else:
                sys.exit("Impossible3! The cigar error: {0}\n{1}".format(read.cigarstring,read.to_string()))
        if not (len(MDrest) == 0 or (len(MDrest) == 1 and MDrest[0] == 0)):
            sys.stderr.write("{0} has the error MD flag".format(read.query_name))
            sys.stderr.flush()
            sys.exit("please check! mark20")
        blocks.append(block(read,qstart,qend,cstart,cend,mycigar,uniqstarsingleread(read),raw_seq,raw_qual,raw_cigar,MDcur,read.get_tag("MD")))
        MDcur = ""
    return blocks

import copy
def printPartBlock(read,blocks,read1_or_read2,fout,pet_index):
    Npart = 0
    if len(blocks) < 2:
        return pet_index,Npart
    for i in range(len(blocks)-1):
        Npart += 1
        pet_index += 1
        this_block = blocks[i]
        next_block = blocks[i+1]
        raw_seq = read.get_forward_sequence()
        raw_qual = pysam.qualities_to_qualitystring(read.get_forward_qualities())
        raw_cigar = read.cigarstring
        pairseq = raw_seq[this_block.start:next_block.end]
        pairqual = raw_qual[this_block.start:next_block.end]
        strand = "+"
        if read.is_reverse:
            strand = "-"
        plus_or_minus = None
        if (strand == "+" and read1_or_read2 == "2") or (strand == "-" and read1_or_read2 == "1"):
            plus_or_minus = "Plus"
        elif (strand == "-" and read1_or_read2 == "2") or (strand == "+" and read1_or_read2 == "1"):
            plus_or_minus = "Minus"
        else:
            sys.exit("strand error")
        
        pair1 = copy.deepcopy(read)
        pair1.query_name = "Part_{0}_{1}_{2}".format(pet_index,plus_or_minus,pair1.query_name)
        pair1.set_tag("RS",raw_seq)
        pair1.set_tag("RQ",raw_qual)
        pair1.set_tag("RC",raw_cigar)
        pair1.cigartuples = this_block.newcigar
        pair1.reference_start = this_block.rstart
        pair1.set_tag("RD",this_block.raw_MD)
        pair1.set_tag("MD",this_block.MD)
        if strand == "+":
            pair1.cigarstring = "{0}{1}S".format(pair1.cigarstring,next_block.end-next_block.start + coutFragmentDis(this_block.start,this_block.end,next_block.start,next_block.end))
            pair1.seq = pairseq
            pair1.qual = pairqual
        else:
            pair1.cigarstring = "{0}S{1}".format(next_block.end-next_block.start + coutFragmentDis(this_block.start,this_block.end,next_block.start,next_block.end),pair1.cigarstring)
            pair1.seq = rc(pairseq)
            pair1.qual = pairqual[::-1]
        
        pair2 = copy.deepcopy(read)
        pair2.query_name = "Part_{0}_{1}_{2}".format(pet_index,plus_or_minus,pair2.query_name)
        pair2.set_tag("RS",raw_seq)
        pair2.set_tag("RQ",raw_qual)
        pair2.set_tag("RC",raw_cigar)
        pair2.cigartuples = next_block.newcigar
        pair2.reference_start = next_block.rstart
        pair2.set_tag("RD",next_block.raw_MD)
        pair2.set_tag("MD",next_block.MD)
        if strand == "+":
            pair2.cigarstring = "{0}S{1}".format(this_block.end-this_block.start + coutFragmentDis(this_block.start,this_block.end,next_block.start,next_block.end),pair2.cigarstring)
            pair2.seq = pairseq
            pair2.qual = pairqual
        else:
            pair2.cigarstring = "{0}{1}S".format(pair2.cigarstring,this_block.end-this_block.start + coutFragmentDis(this_block.start,this_block.end,next_block.start,next_block.end))
            pair2.seq = rc(pairseq)
            pair2.qual = pairqual[::-1]
        if strand == "+":
            fout.write(pair1)
            fout.write(pair2)
        else:
            fout.write(pair2)
            fout.write(pair1)
    return pet_index,Npart

def printChimericBlock(reads,readsblocks,read1_or_read2,fout,pet_index):
    Nchimeric = 1
    pet_index += 1

    raw_seq = reads[0].get_forward_sequence()
    raw_qual = pysam.qualities_to_qualitystring(reads[0].get_forward_qualities())

    strand1 = "+"
    if reads[0].is_reverse:
        strand1 = "-"
    if (strand1 == "+" and read1_or_read2 == "2") or (strand1 == "-" and read1_or_read2 == "1"):
        plus_or_minus1 = "Plus"
    elif (strand1 == "-" and read1_or_read2 == "2") or (strand1 == "+" and read1_or_read2 == "1"):
        plus_or_minus1 = "Minus"
    else:
        sys.exit("strand error")
    pair1 = copy.deepcopy(reads[0])
    if len(readsblocks[0]) == 1:
        pair1.query_name = "ChimericWhole_{0}_{1}_Head_{2}".format(pet_index,plus_or_minus1,pair1.query_name)
    else:
        pair1.query_name = "ChimericPart_{0}_{1}_Head_{2}".format(pet_index,plus_or_minus1,pair1.query_name)
    pair1.set_tag("RS",raw_seq)
    pair1.set_tag("RQ",raw_qual)
    raw_cigar1 = reads[0].cigarstring
    pair1.set_tag("RC",raw_cigar1)
    pair1.cigartuples = readsblocks[0][-1].newcigar
    pair1.reference_start = readsblocks[0][-1].rstart
    this_block = readsblocks[0][-1]
    pair1.set_tag("RD",this_block.raw_MD)
    pair1.set_tag("MD",this_block.MD)
    
    strand2 = "+"
    if reads[1].is_reverse:
        strand2 = "-"
    if (strand2 == "+" and read1_or_read2 == "2") or (strand2 == "-" and read1_or_read2 == "1"):
        plus_or_minus2 = "Plus"
    elif (strand2 == "-" and read1_or_read2 == "2") or (strand2 == "+" and read1_or_read2 == "1"):
        plus_or_minus2 = "Minus"
    else:
        sys.exit("strand error")
    pair2 = copy.deepcopy(reads[1])
    if len(readsblocks[1]) == 1:
        pair2.query_name = "ChimericWhole_{0}_{1}_Tail_{2}".format(pet_index,plus_or_minus2,pair2.query_name)
    else:
        pair2.query_name = "ChimericPart_{0}_{1}_Tail_{2}".format(pet_index,plus_or_minus2,pair2.query_name)
    pair2.set_tag("RS",raw_seq)
    pair2.set_tag("RQ",raw_qual)
    raw_cigar2 = reads[1].cigarstring
    pair2.set_tag("RC",raw_cigar2)
    pair2.cigartuples = readsblocks[1][0].newcigar
    pair2.reference_start = readsblocks[1][0].rstart
    next_block = readsblocks[1][0]
    pair2.set_tag("RD",next_block.raw_MD)
    pair2.set_tag("MD",next_block.MD)

    pairseq = raw_seq[this_block.start:next_block.end]
    pairqual = raw_qual[this_block.start:next_block.end]

    if strand1 == "+":
        pair1.cigarstring = "{0}{1}S".format(pair1.cigarstring,next_block.end-next_block.start + coutFragmentDis(this_block.start,this_block.end,next_block.start,next_block.end))
        pair1.seq = pairseq
        pair1.qual = pairqual
    else:
        pair1.cigarstring = "{0}S{1}".format(next_block.end-next_block.start + coutFragmentDis(this_block.start,this_block.end,next_block.start,next_block.end),pair1.cigarstring)
        pair1.seq = rc(pairseq)
        pair1.qual = pairqual[::-1]
    if strand2 == "+":
        pair2.cigarstring = "{0}S{1}".format(this_block.end-this_block.start + coutFragmentDis(this_block.start,this_block.end,next_block.start,next_block.end),pair2.cigarstring)
        pair2.seq = pairseq
        pair2.qual = pairqual
    else:
        pair2.cigarstring = "{0}{1}S".format(pair2.cigarstring,this_block.end-this_block.start + coutFragmentDis(this_block.start,this_block.end,next_block.start,next_block.end))
        pair2.seq = rc(pairseq)
        pair2.qual = pairqual[::-1]
    fout.write(pair1)
    fout.write(pair2)
    return pet_index, Nchimeric

def getPairFromChimeric(reads,read1_or_read2,fout,pet_index):
    readsblocks = []
    Npart = 0
    Nchimeric = 0
    # print Part
    for read in reads:
        blocks = starRead2Blocks(read)
        readsblocks.append(blocks)
        pet_index,Ntmp = printPartBlock(read,blocks,read1_or_read2,fout,pet_index)
        Npart += Ntmp
    # print chimeric
    if len(readsblocks) != 2 or len(reads) != 2:
        for ii in reads:
            print(ii)
        sys.exit("getPairFromChimeric error!")
    pet_index,Ntmp = printChimericBlock(reads,readsblocks,read1_or_read2,fout,pet_index)
    Nchimeric += Ntmp
    return pet_index, Npart, Nchimeric

def getPairFromStarReads(reads,read1_or_read2,fout,pet_index):
    readsblocks = []
    Npart = 0
    Nchimeric = 0
    # print Part
    for read in reads:
        blocks = starRead2Blocks(read)
        readsblocks.append(blocks)
        pet_index,Ntmp = printPartBlock(read,blocks,read1_or_read2,fout,pet_index)
        Npart += Ntmp
    return pet_index, Npart, Nchimeric

def bwaRead2Blocks(reads):
    blocks = []
    if reads[0].is_unmapped:
        return blocks
    raw_seq = None
    raw_qual = None
    for ii in reads:
        if not ii.is_supplementary:
            raw_seq = ii.get_forward_sequence()
            raw_qual = pysam.qualities_to_qualitystring(ii.get_forward_qualities())
            break
    for ii in reads:
        s, e, newcigar = getReadPosPro(ii)
        blocks.append(block(ii,s,e,ii.reference_start,ii.reference_end,newcigar,mapstatbwa(ii),raw_seq,raw_qual,ii.cigarstring,ii.get_tag("MD"),ii.get_tag("MD")))
    blocks.sort(key=lambda x:(x.start,x.end))
    return blocks

def fltBwaBlocks(blocks,nfilter=2):
    fltblocks = []
    for ii in blocks:
        if ii.mapstat == "uniqmap":
            fltblocks.append(ii)
    if len(fltblocks) < nfilter:
        return "uniqmap_block_less_than_{0}".format(nfilter),None

    blocks = []
    for i in range(len(fltblocks)):
        conflict_i = 0
        for j in range(len(fltblocks)):
            if i == j:
                continue
            ccb = coutFragmentDis(fltblocks[i].start,fltblocks[i].end,fltblocks[j].start,fltblocks[j].end)
            if ccb < 0:
                conflict_i -= ccb
        if conflict_i < 0.25*(fltblocks[i].end - fltblocks[i].start): # 重叠的reads数目太多的也不要
            blocks.append(fltblocks[i])
    if len(blocks) < nfilter:
        return "overlap_more",None

    effective_mapped = set()
    for i in range(len(blocks)):
        for j in range(blocks[i].start,blocks[i].end):
            effective_mapped.add(j)
    effective_len = len(effective_mapped)
    #### 这里是有潜在bug的。如果reads断成3端，其中有一段因为某种原因被过滤掉了，那么片段的加和就很难超越80%了。
    #### 所以正确的做法应该是，每次临近的两个片段来做判断。 好的：------ --------  坏的： -----                   -----                         -------
    if effective_len < 0.8*len(blocks[0].raw_seq): # 不管是前面的0.25还是这里的0.8都是经验值
        return "gap_is_large",None
    
    if len(blocks) > 4:
        return "remove_over_{0}_blocks".format(len(blocks)),None
    return "pass",blocks

def printBwaBlock(blocks,read1_or_read2,fout,pet_index):
    Nchimeric = 0
    Npart = 0
    if len(blocks) < 2:
        return pet_index, Npart, Nchimeric
    for i in range(len(blocks)-1):
        this_read = blocks[i].read
        this_strand = "+"
        if this_read.is_reverse:
            this_strand = "-"
        if (this_strand == "+" and read1_or_read2 == "2") or (this_strand == "-" and read1_or_read2 == "1"):
            plus_or_minus_this = "Plus"
        elif (this_strand == "-" and read1_or_read2 == "2") or (this_strand == "+" and read1_or_read2 == "1"):
            plus_or_minus_this = "Minus"
        else:
            sys.exit("chrck2: {0}\t{1}".format(blocks[i].to_string(),read1_or_read2))
            
        next_read = blocks[i+1].read # 在这里，只考虑了连续的两块，不考虑间接的两块
        next_strand = "+"
        if next_read.is_reverse:
            next_strand = "-"
        if (next_strand == "+" and read1_or_read2 == "2") or (next_strand == "-" and read1_or_read2 == "1"):
                plus_or_minus_next = "Plus"
        elif (next_strand == "-" and read1_or_read2 == "2") or (next_strand == "+" and read1_or_read2 == "1"):
            plus_or_minus_next = "Minus"
        else:
            sys.exit("chrck3: {0}\t{1}".format(next_read.to_string(),read1_or_read2))
        
        # 判断是否越界其实是可以不做的。没有太大的必要。因为reads可能和参考有一定的差异。同时现在不能判断延申方向
        
        is_singleton = False
        if this_strand == next_strand and this_read.reference_name == next_read.reference_name:
            if this_strand == "+" and next_read.reference_start > this_read.reference_end:	#singleton
                is_singleton = True
            elif this_strand == "-" and next_read.reference_end < this_read.reference_start:	#singleton
                is_singleton = True

        pet_index += 1
        pair1 = copy.deepcopy(this_read)
        pair1.set_tag("RS",blocks[i].raw_seq)
        pair1.set_tag("RQ",blocks[i].raw_qual)
        pair1.set_tag("RG",this_read.cigarstring)
        pair2 = copy.deepcopy(next_read)
        pair2.set_tag("RS",blocks[i+1].raw_seq)
        pair2.set_tag("RQ",blocks[i+1].raw_qual)
        pair2.set_tag("RG",next_read.cigarstring)

        pair1.set_tag("RD",blocks[i].raw_MD)
        pair1.set_tag("MD",blocks[i].MD)
        pair2.set_tag("RD",blocks[i+1].raw_MD)
        pair2.set_tag("MD",blocks[i+1].MD)

        pair1.seq = blocks[i].raw_seq[blocks[i].start:blocks[i+1].end]
        pair1.qual = blocks[i].raw_qual[blocks[i].start:blocks[i+1].end]
        pair1.cigartuples = blocks[i].newcigar
        pair1.reference_start = blocks[i].rstart
        pair2.seq = pair1.seq
        pair2.qual = pair1.qual
        pair2.cigartuples = blocks[i+1].newcigar
        pair2.reference_start = blocks[i+1].rstart
        if is_singleton:
            Npart += 1
            pair1.query_name = "Part_{0}_{1}_{2}".format(pet_index,plus_or_minus_this,pair1.query_name)
            pair2.query_name = "Part_{0}_{1}_{2}".format(pet_index,plus_or_minus_next,pair2.query_name)
            if this_strand == "+" and next_strand == "+":
                pair1.cigarstring = "{0}{1}S".format(pair1.cigarstring,blocks[i+1].end - blocks[i+1].start + coutFragmentDis(blocks[i].start,blocks[i].end,blocks[i+1].start,blocks[i+1].end))
                pair2.cigarstring = "{0}S{1}".format(blocks[i].end - blocks[i].start + coutFragmentDis(blocks[i].start,blocks[i].end,blocks[i+1].start,blocks[i+1].end),pair2.cigarstring)
                fout.write(pair1)
                fout.write(pair2)
            elif this_strand == "-" and next_strand == "-":
                pair1.cigarstring = "{0}S{1}".format(blocks[i+1].end - blocks[i+1].start + coutFragmentDis(blocks[i].start,blocks[i].end,blocks[i+1].start,blocks[i+1].end),pair1.cigarstring)
                pair2.cigarstring = "{0}{1}S".format(pair2.cigarstring,blocks[i].end - blocks[i].start + coutFragmentDis(blocks[i].start,blocks[i].end,blocks[i+1].start,blocks[i+1].end))
                tmp = pair1.qual[::-1]
                pair1.seq = rc(pair1.seq)
                pair1.qual = tmp
                tmp = pair2.qual[::-1]
                pair2.seq = rc(pair2.seq)
                pair2.qual = tmp
                fout.write(pair2)
                fout.write(pair1)
            else:
                sys.exit("error1! \n{0}\n{1}\n".format(this_read.to_string(),next_read.to_string()))
        else:
            Nchimeric += 1
            if this_strand == "+":
                pair1.cigarstring = "{0}{1}S".format(pair1.cigarstring,blocks[i+1].end - blocks[i+1].start + coutFragmentDis(blocks[i].start,blocks[i].end,blocks[i+1].start,blocks[i+1].end))
            else:
                tmp = pair1.qual[::-1]
                pair1.seq = rc(pair1.seq)
                pair1.qual = tmp
                pair1.cigarstring = "{0}S{1}".format(blocks[i+1].end - blocks[i+1].start + coutFragmentDis(blocks[i].start,blocks[i].end,blocks[i+1].start,blocks[i+1].end),pair1.cigarstring)
            if next_strand == "+":
                pair2.cigarstring = "{0}S{1}".format(blocks[i].end - blocks[i].start + coutFragmentDis(blocks[i].start,blocks[i].end,blocks[i+1].start,blocks[i+1].end), pair2.cigarstring)
            else:
                tmp = pair2.qual[::-1]
                pair2.seq = rc(pair2.seq)
                pair2.qual = tmp
                pair2.cigarstring = "{0}{1}S".format(pair2.cigarstring,blocks[i].end - blocks[i].start + coutFragmentDis(blocks[i].start,blocks[i].end,blocks[i+1].start,blocks[i+1].end))
            pair1.query_name = "ChimericPart_{0}_{1}_Head_{2}".format(pet_index,plus_or_minus_this,pair1.query_name)
            fout.write(pair1)
            pair2.query_name = "ChimericPart_{0}_{1}_Tail_{2}".format(pet_index,plus_or_minus_next,pair2.query_name)
            fout.write(pair2)
    return pet_index, Npart, Nchimeric

def readChimeric(file,taginfo):
    fchim = SamRead(file)
    star_chimeric_full = {}
    star_chimeric_notfull = {}
    for reads in fchim:
        addGroup(reads,"GI",taginfo)
        if len(reads) != 2:
            for ii in reads:
                print(ii)
            sys.exit("error! The STAR chimeric file's reads has more/less than 2 split fragment! \n")
        readlen,beforeS,afterS,maplen,dis = checkChimeric(reads)
        if beforeS < 15 and afterS < 15:
            star_chimeric_full[reads[0].query_name] = reads
        else:
            star_chimeric_notfull[reads[0].query_name] = reads
    sys.stdout.write("full chimeric: {0}\n".format(len(star_chimeric_full)))
    sys.stdout.write("not full chimeric: {0}\n".format(len(star_chimeric_notfull)))
    sys.stdout.flush()
    return star_chimeric_full,star_chimeric_notfull

def getPairFromBwaReads(reads,read1_or_read2,fout,pet_index):
    Npart = 0
    Nchimeric = 0
    flag, blocks = fltBwaBlocks(bwaRead2Blocks(reads))
    # print Part
    if flag != "pass":
        return flag,pet_index,0,0
    pet_index,Npart,Nchimeric = printBwaBlock(blocks,read1_or_read2,fout,pet_index)
    return flag,pet_index, Npart, Nchimeric

def detstarbwa(star_chimeric_full,star_chimeric_notfull,starreads,bwareads,mapflag):
    mapflag.append("##")
    tmpstarflag = checkStarAlign(starreads)
    mapflag.append(tmpstarflag[0])
    tmpbwaflag = checkBwaAlign(bwareads)
    mapflag.append(tmpbwaflag[0])
    # print(tmpstarflag)
    # print(tmpbwaflag)

    mapflag.extend(["##","##"])
    mapflag.extend(checkSB(starreads,bwareads))
    if starreads[0].query_name in star_chimeric_full:
        mapflag[1] = "chimeric#full"
        checkresult = checkCSB(star_chimeric_full[starreads[0].query_name],starreads,bwareads)
        mapflag[4] = checkresult[0]
        mapflag[5] = checkresult[1]
        return star_chimeric_full[starreads[0].query_name],"chimeric","3"
        # if "multi" in mapflag[2]: # STAR比对的chimeric的结果只会报道两个比对位置，不会报道uniq还是multi的。所以需要用STAR的align的结果来判断。
        #     if "unmap#full" == mapflag[3]:
        #         return None,None,"1"
        #     else:
        #         # return bwareads,"bwa","2"
        #         return star_chimeric_full[starreads[0].query_name],"chimeric","2"
        # else:
        #     return star_chimeric_full[starreads[0].query_name],"chimeric","3"
    elif starreads[0].query_name in star_chimeric_notfull:
        mapflag[1] = "chimeric#notfull"
        checkresult = checkCSB(star_chimeric_notfull[starreads[0].query_name],starreads,bwareads)
        mapflag[4] = checkresult[0]
        mapflag[5] = checkresult[1]
        return star_chimeric_notfull[starreads[0].query_name],"chimeric","7"
        # if "multi" in mapflag[2]: # STAR比对的chimeric的结果只会报道两个比对位置，不会报道uniq还是multi的。所以需要用STAR的align的结果来判断。
        #     if "unmap#full" == mapflag[3]:
        #         return None,None,"4"
        #     else:
        #         # return bwareads,"bwa","5"
        #         return star_chimeric_notfull[starreads[0].query_name],"chimeric","5"
        # else:
        #     if "unmap#full" == mapflag[3]:
        #         return star_chimeric_notfull[starreads[0].query_name],"chimeric","6"
        #     else:
        #         # return bwareads,"bwa","7"
        #         return star_chimeric_notfull[starreads[0].query_name],"chimeric","7"
    else:
        mapflag[1] = "nochimeric"
        if "uniqmap#full" == mapflag[2]:
            return starreads,"star","8"
        elif "multimap#full" == mapflag[2]:
            return None,None,"9"
        elif "uniqmap#notfull" == mapflag[2]:
            if tmpbwaflag[2] + tmpbwaflag[1] < tmpstarflag[2] + tmpstarflag[1] + 15: # 如果bwa拥有比star长15bp以上的比对结果，选择bwa的比对结果来考虑chimeric
                return bwareads,"bwa","10"
            else:
                return starreads,"star","11"
        elif "multimap#notfull" == mapflag[2]:
            if tmpbwaflag[2] + tmpbwaflag[1] < tmpstarflag[2] + tmpstarflag[1] + 15:
                return bwareads,"bwa","12"
            else:
                return None,None,"13"
        else:
            if "unmap#full" == mapflag[3]:
                return None,None,"14"
            else:
                return bwareads,"bwa","15"

def detstar(star_chimeric_full,star_chimeric_notfull,starreads,mapflag):
    mapflag.append("##")
    tmpstarflag = checkStarAlign(starreads)
    mapflag.append(tmpstarflag[0])
    mapflag.append("-")

    mapflag.extend(["##","##"])
    mapflag.extend(["##","##"])
    if starreads[0].query_name in star_chimeric_full:
        mapflag[1] = "chimeric#full"
        checkresult = checkCS(star_chimeric_full[starreads[0].query_name],starreads)
        mapflag[4] = checkresult[0]
        mapflag[5] = checkresult[1]
        return star_chimeric_full[starreads[0].query_name],"chimeric","3"
        # if "multi" in mapflag[2]: # STAR比对的chimeric的结果只会报道两个比对位置，不会报道uniq还是multi的。所以需要用STAR的align的结果来判断。
        #     if "unmap#full" == mapflag[3]:
        #         return star_chimeric_full[starreads[0].query_name],"chimeric","1"
        #     else:
        #         # return bwareads,"bwa","2"
        #         return star_chimeric_full[starreads[0].query_name],"chimeric","2"
        # else:
        #     return star_chimeric_full[starreads[0].query_name],"chimeric","3"
    elif starreads[0].query_name in star_chimeric_notfull:
        mapflag[1] = "chimeric#notfull"
        checkresult = checkCS(star_chimeric_notfull[starreads[0].query_name],starreads)
        mapflag[4] = checkresult[0]
        mapflag[5] = checkresult[1]
        return star_chimeric_notfull[starreads[0].query_name],"chimeric","7"
        # if "multi" in mapflag[2]: # STAR比对的chimeric的结果只会报道两个比对位置，不会报道uniq还是multi的。所以需要用STAR的align的结果来判断。
        #     if "unmap#full" == mapflag[3]:
        #         # return None,None,"4"
        #         return star_chimeric_notfull[starreads[0].query_name],"chimeric","4"
        #     else:
        #         # return bwareads,"bwa","5"
        #         return star_chimeric_notfull[starreads[0].query_name],"chimeric","5"
        # else:
        #     if "unmap#full" == mapflag[3]:
        #         return star_chimeric_notfull[starreads[0].query_name],"chimeric","6"
        #     else:
        #         # return bwareads,"bwa","7"
        #         return star_chimeric_notfull[starreads[0].query_name],"chimeric","7"
    else:
        mapflag[1] = "nochimeric"
        if "uniqmap#full" == mapflag[2]:
            return starreads,"star","8"
        elif "multimap#full" == mapflag[2]:
            return None,None,"9"
        elif "uniqmap#notfull" == mapflag[2]:
            return starreads,"star","11"
        elif "multimap#notfull" == mapflag[2]:
            return None,None,"13"
        else:
            return None,None,"14"

def detsamechimeric(ff,rr):
    R1blocks = []
    for read in ff:
        R1blocks.extend(starRead2Blocks(read))
    R2blocks = []
    for read in rr:
        R2blocks.extend(starRead2Blocks(read))
    sameset = set()
    for block1 in R1blocks:
        for block2 in R2blocks:
            if block1.read.reference_name == block2.read.reference_name:
                if block1.rstart - block2.rstart > -100 and block1.rstart - block2.rstart < 100:
                    sameset.add(block1)
                    break
                if block1.rend - block2.rend > -100 and block1.rend - block2.rend < 100:
                    sameset.add(block1)
                    break
    if len(sameset) > len(R1blocks)*0.8:
        if len(R1blocks[0].raw_seq) >= len(R2blocks[0].raw_seq):
            return "same","R1"
        else:
            return "same","R2"
    elif len(sameset) > len(R1blocks)*0.4:
        if len(R1blocks[0].raw_seq) >= len(R2blocks[0].raw_seq):
            return "similar","R1"
        else:
            return "similar","R2"
    else:
        return "diff","R1R2"

def detsamebwa(ff,rr):
    R1blocks = fltBwaBlocks(bwaRead2Blocks(ff))
    R2blocks = fltBwaBlocks(bwaRead2Blocks(rr))
    sameset = set()
    for block1 in R1blocks:
        for block2 in R2blocks:
            if block1.read.reference_name == block2.read.reference_name:
                if block1.rstart - block2.rstart > -100 and block1.rstart - block2.rstart < 100:
                    sameset.add(block1)
                    break
                if block1.rend - block2.rend > -100 and block1.rend - block2.rend < 100:
                    sameset.add(block1)
                    break
    if len(sameset) > len(R1blocks)*0.8:
        if len(R1blocks[0].raw_seq) >= len(R2blocks[0].raw_seq):
            return "same","R1"
        else:
            return "same","R2"
    elif len(sameset) > len(R1blocks)*0.4:
        if len(R1blocks[0].raw_seq) >= len(R2blocks[0].raw_seq):
            return "similar","R1"
        else:
            return "similar","R2"
    else:
        return "diff","R1R2"

def checkReadsID(starreads,bwareads):
    if starreads == None or bwareads == None:
        sys.stderr.write("The reads count of star align file and bwa algin file don't match\n")
        sys.stderr.flush()
        sys.exit("{0}\t{1}\n".format(starreads[0].query_name,bwareads[0].query_name))
    elif starreads[0].query_name != bwareads[0].query_name:
        sys.stderr.write("The reads of star align file and bwa algin file don't match\n")
        sys.stderr.flush()
        sys.exit("{0}\t{1}\n".format(starreads[0].query_name,bwareads[0].query_name))

def updateDict3(mystatic,Npart,Nchimeric, Npair):
    if "Part" not in mystatic:
        mystatic["Part"] = {}
    if Npart not in mystatic["Part"]:
        mystatic["Part"][Npart] = 0
    mystatic["Part"][Npart] += 1
    if "Chimeric" not in mystatic:
        mystatic["Chimeric"] = {}
    if Nchimeric not in mystatic["Chimeric"]:
        mystatic["Chimeric"][Nchimeric] = 0
    mystatic["Chimeric"][Nchimeric] += 1
    if "PairAlign" not in mystatic:
        mystatic["PairAlign"] = {}
    if Npair not in mystatic["PairAlign"]:
        mystatic["PairAlign"][Npair] = 0
    mystatic["PairAlign"][Npair] += 1

def updateDict(mystatic,Npart,Nchimeric):
    if "Part" not in mystatic:
        mystatic["Part"] = {}
    if Npart not in mystatic["Part"]:
        mystatic["Part"][Npart] = 0
    mystatic["Part"][Npart] += 1
    if "Chimeric" not in mystatic:
        mystatic["Chimeric"] = {}
    if Nchimeric not in mystatic["Chimeric"]:
        mystatic["Chimeric"][Nchimeric] = 0
    mystatic["Chimeric"][Nchimeric] += 1

def updateDict2(mystatic,flag):
    if flag not in mystatic:
        mystatic[flag] = 0
    mystatic[flag] += 1

class SamRead:
    """
    get sam read
    """
    def __init__(self,infile):
        self.fin = readSam(infile)
    def __iter__(self):
        flag = None
        outlist = []
        for read in self.fin:
            if flag != None and flag != read.query_name:
                yield outlist
                outlist = []
            outlist.append(read)
            flag = read.query_name
        self.fin.close()
        yield outlist
    def header(self):
        return self.fin.header

class FaRead:
    def __init__(self,infile):
        if infile == None or infile == "-":
            self.fin = sys.stdin
        else:
            self.fin = readFile(infile)
    def __iter__(self):
        flag = None
        outstr = ""
        for line in self.fin:
            if flag == None:
                flag = line.strip()[1:]
            elif line.startswith(">"):
                yield flag,outstr
                flag = line.strip()[1:]
                outstr = ""
            else:
                outstr += line.strip()
        self.fin.close()
        yield flag,outstr
        
class FqRead:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        for line in self.fin:
            rid = line.strip().split()[0]
            rid = rid[1:]
            rseq = self.fin.readline()
            rseq = rseq.strip()
            rsyb = self.fin.readline()
            rsyb = rsyb.strip()
            rqual = self.fin.readline()
            rqual = rqual.strip()
            yield rid,rseq,rsyb,rqual
        self.fin.close()

class FqRead2:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        for line in self.fin:
            rid = line.strip()
            rseq = self.fin.readline()
            rseq = rseq.strip()
            rsyb = self.fin.readline()
            rsyb = rsyb.strip()
            rqual = self.fin.readline()
            rqual = rqual.strip()
            yield rid,rseq,rsyb,rqual
        self.fin.close()

class readseq:
    def __init__(self,rid,seq,length,qual,strand):
        self.rid = rid
        self.seq = seq
        self.length = length
        self.qual = qual
        self.strand = strand
        self.readinfo = []
    def check(self):
        if len(self.seq) != self.length:
            sys.exit("you read have error{0}".format(self.__str__()))
        if len(self.qual) != self.length:
            sys.exit("you read have error{0}".format(self.__str__()))
    def __str__(self):
        return "{0}\t{1}\t{2}".format(self.seq,self.qual,self.qual)
    def __repr__(self):
        return self.__str__()

class readLoc:
    def __init__(self,read,start,end):
        self.read = read
        self.start = start
        self.end = end
    def __str__(self):
        return "{0}\t{1}\t{2}".format(self.start,self.end,self.read.to_string())
    def __repr__(self):
        return self.__str__()

def LocalAlign(aa,bb,start = 0):
    MismatchScore = -1
    MatchScore = 1
    IndelScore = -2 # 设置的更大，就可以避免通过indel
    maxI = 0
    maxJ = 0
    minI = 0
    minJ = 0
    m = len(aa)
    n = len(bb)
    scoreMatrix = []
    maxScore = 0
    for i in range(m+1):
        scoreMatrix.append([])
        for j in range(n+1):
            scoreMatrix[i].append(0)
        
    for i in range(m):
        for j in range(start,n):
            if aa[i] == bb[j]:
                scoreMatrix[i+1][j+1] = scoreMatrix[i][j]+MatchScore
            else:
                scoreMatrix[i+1][j+1] = scoreMatrix[i][j]+MismatchScore
            if scoreMatrix[i+1][j+1] < scoreMatrix[i+1][j] + IndelScore:
                scoreMatrix[i+1][j+1] = scoreMatrix[i+1][j]+IndelScore
            if scoreMatrix[i+1][j+1] < scoreMatrix[i][j+1] + IndelScore:
                scoreMatrix[i+1][j+1] = scoreMatrix[i][j+1]+IndelScore
            if scoreMatrix[i+1][j+1] < 0:
                scoreMatrix[i+1][j+1] = 0
            if maxScore < scoreMatrix[i+1][j+1]:
                maxScore = scoreMatrix[i+1][j+1]
                maxI = i+1
                maxJ = j+1
    # for i in scoreMatrix:
    #     print(i)
    # print(maxScore)
        
    i = maxI-1
    j = maxJ-1
    nMatches = 0
    nMismatches = 0
    nIndels = 0
    alignedStatus = ""
    str1 = ""
    str2 = ""
    scoreGreaterThan0 = 1

    while m-1 > i:
        str1 += aa[m-1]
        str2 += "-"
        alignedStatus += " "
        m -= 1
    while n-1 > j:
        str1 += "-"
        str2 += bb[n-1]
        alignedStatus += " "
        n -= 1
    while i>=0 or j>=0:
        if i < 0:
            str1 += "-"
            str2 += bb[j]
            alignedStatus += " "
            j -= 1
            if scoreGreaterThan0 == 1:
                nIndels += 1
        elif j < 0:
            str1 += aa[i]
            str2 += "-"
            alignedStatus += " "
            i -= 1
            if scoreGreaterThan0 == 1:
                nIndels += 1
        elif scoreMatrix[i+1][j+1] == scoreMatrix[i][j] + MatchScore and aa[i] == bb[j]:
            str1 += aa[i]
            str2 += bb[j]
            alignedStatus += "|"
            i -= 1
            j -= 1
            if scoreGreaterThan0 == 1:
                nMatches += 1
        elif scoreMatrix[i+1][j+1] == scoreMatrix[i][j] + MismatchScore:
            str1 += aa[i]
            str2 += bb[j]
            alignedStatus += "X"
            i -= 1
            j -= 1
            if scoreGreaterThan0 == 1:
                nMismatches += 1
        elif scoreMatrix[i+1][j+1] == scoreMatrix[i+1][j] + IndelScore:
            str1 += "-"
            str2 += bb[j]
            alignedStatus += " "
            j -= 1
            if scoreGreaterThan0 == 1:
                nIndels += 1
        else:
            str1 += aa[i]
            str2 += "-"
            alignedStatus += " "
            i -= 1
            if scoreGreaterThan0 == 1:
                nIndels += 1
        if scoreGreaterThan0 == 1:
            if scoreMatrix[i+1][j+1] <= 0:
                minI = i+1
                minJ = j+1
                scoreGreaterThan0 = 0
    return maxScore,nMatches,nMismatches,nIndels,minI,maxI,minJ,maxJ,str1[::-1],alignedStatus[::-1],str2[::-1]

from bx.intervals.intersection import Intersecter, Interval
def regionTree(resFrag,chrom,start,end,flag=None):
    if chrom not in resFrag:
        resFrag[chrom] = Intersecter()
    resFrag[chrom].add_interval(Interval(int(start),int(end),flag))

def regionFind(tree,chrom,start,end):
    if chrom not in tree:
        return []
    return tree[chrom].find(start,end)

def checkPairtagRead(read1,read2):
    tmp1 = read1.query_name.split("_")
    tmp2 = read2.query_name.split("_")
    if tmp1[1] != tmp2[1]:
        return True
    if tmp1[-1] != tmp2[-1]:
        return True
    return False

def getStrand(read):
    info = read.query_name.split("_")
    if info[2] == "Plus":
        return "+"
    elif info[2] == "Minus":
        return "-"
    else:
        sys.stderr.write("strand error {0}".format(read.tostring()))
        sys.exit(1)

from collections import OrderedDict
def parseStringgtf(mystr):
    d = OrderedDict()
    if mystr.endswith(";"):
        mystr = mystr[:-1]
    for i in mystr.strip().split(";"):
        tmp = i.strip().split()
        if tmp == "":
            continue
        d[tmp[0]] = tmp[1].replace("\"","")
    return d
