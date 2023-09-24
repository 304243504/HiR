import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import readFile,regionTree,readSam,checkPairtagRead,regionFind,writeSam
import argparse

"""
Op  BAM Description Consumes_query  Consumes_reference
M 0 alignment match (can be a sequence match or mismatch) yes yes
I 1 insertion to the reference yes no
D 2 deletion from the reference no yes
N 3 skipped region from the reference no yes
S 4 soft clipping (clipped sequences present in SEQ) yes no
H 5 hard clipping (clipped sequences NOT present in SEQ) no no
P 6 padding (silent deletion from padded reference) no no
= 7 sequence match yes yes
X 8 sequence mismatch yes yes
"""

parser = argparse.ArgumentParser(description = "split interaction file into intra-gene, inter-gene, without-gene interaction")
parser.add_argument("-g","--generegion",required=True,help="the input gene region file")
parser.add_argument("-i","--input",required=True,help="the input interaction files, comma separated (file1,file2,file3,...  )")
parser.add_argument("-o","--outputprefix",required=True,help="the output prefix")
args = parser.parse_args()

if __name__ == "__main__":
    it = {} # Interval Tree
    fin = readFile(args.generegion)
    for line in fin:
        tmp = line.strip().split("\t")
        if tmp[4] == "gene":
            regionTree(it,tmp[0],tmp[1],tmp[2],tmp[3:])
    fin.close()
    
    files = args.input.split(",")
    fin = readSam(files[0])
    fintra = writeSam("{0}.interaction.intraMolecular.bam".format(args.outputprefix),fin.header)
    finter = writeSam("{0}.interaction.interMolecular.bam".format(args.outputprefix),fin.header)
    fin.close()
    nintra = 0
    nwith = 0
    nwithout = 0
    nwithhalf = 0
    myhash = {}

    for tmpfile in files:
        fin = readSam(tmpfile)
        for read1 in fin:
            read2 = next(fin)
            if read1.is_reverse:
                strand1 = "-"
            else:
                strand1 = "+"
            if read2.is_reverse:
                strand2 = "-"
            else:
                strand2 = "+"
            hit1 = regionFind(it,read1.reference_name,read1.reference_start,read1.reference_end)
            hit2 = regionFind(it,read2.reference_name,read2.reference_start,read2.reference_end)
            geneset1 = set()
            geneset2 = set()
            for ii in hit1:
                tmp1 = ii.value
                if strand1 == tmp1[2]:
                    geneset1.add(tmp1[0])
            for jj in hit2:
                tmp2 = jj.value
                if strand2 == tmp2[2]:
                    geneset2.add(tmp2[0])
            read1.set_tag("HG",";".join(sorted(geneset1))) # hit gene
            read2.set_tag("HG",";".join(sorted(geneset2)))
            if len(geneset1.intersection(geneset2)) > 0: # same gene
                read1.set_tag("FG","intra") # flag gene
                read2.set_tag("FG","intra")
                nintra += 1
                fintra.write(read1)
                fintra.write(read2)
            else:
                if len(geneset1) == 0 and len(geneset2) == 0:
                    read1.set_tag("FG","without")
                    read2.set_tag("FG","without")
                    nwithout += 1
                elif len(geneset1) > 0 and len(geneset2) > 0:
                    read1.set_tag("FG","with")
                    read2.set_tag("FG","with")
                    nwith += 1
                else:
                    read1.set_tag("FG","withhalf")
                    read2.set_tag("FG","withhalf")
                    nwithhalf += 1
                finter.write(read1)
                finter.write(read2)

        fin.close()
    fintra.close()
    finter.close()

    print("intra pair tag: {0}".format(nintra))
    print("inter pair tag: {0}".format(nwith+nwithout+nwithhalf))
    print("inter with gene pair tag: {0}".format(nwith))
    print("inter withhalf gene pair tag: {0}".format(nwithhalf))
    print("inter without gene pair tag: {0}".format(nwithout))
    print("--------------------------------------------------------------")


