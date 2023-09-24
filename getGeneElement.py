import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import readFile,writeFile,parseStringgtf
import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser(description = "get gene element")
parser.add_argument("-i","--input",required=True,help="the input file")
parser.add_argument("-o","--output",required=True,help="the output file")
args = parser.parse_args()

class gene:
    def __init__(self,gid,chrom,start,end,strand):
        self.id = gid
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.rna = OrderedDict()
    
    def check(self):
        for ii in self.rna.keys():
            if self.strand != self.rna[ii].strand:
                sys.exit("the gene and transcript has different strand: \n{0}\n{1}\n".format(self.id,ii))

class trans:
    def __init__(self,rid,chrom,start,end,strand):
        self.id = rid
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.exons = []
        self.cds = []
        self.introns = []
        self.cdsintrons = []
    
    def getIntron(self):
        if len(self.exons) > 1:
            self.exons.sort(key=lambda x:(x[0],x[1]))
            for i in range(len(self.exons)-1):
                self.introns.append([self.exons[i][1],self.exons[i+1][0]])
        if len(self.cds) > 1:
            self.cds.sort(key=lambda x:(x[0],x[1]))
            for i in range(len(self.cds)-1):
                self.cdsintrons.append([self.cds[i][1],self.cds[i+1][0]])
    
    def check(self):
        if len(self.cds) > len(self.exons):
            sys.exit("the number of cds is great than exons in {0}".format(self.id))
        if len(self.exons) == 0:
            sys.exit("transcript {0} has no exon\n".format(self.id))
        if self.exons[0][0] != self.start:
            sys.stderr.write("transcript start don't match with exon start {0}\n".format(self.id))
            self.exons[0][0] = self.start
        if self.exons[-1][-1] != self.end:
            sys.stderr.write("transcript end don't match with exon end {0}\n".format(self.id))
            self.exons[-1][-1] = self.end
        if len(self.exons) == 0:
            sys.exit("transcript has no exon {0}".format(self.id))
    
    def getForword(self):
        if self.strand == "-":
            self.exons = self.exons[::-1]
            self.cds = self.cds[::-1]
            self.introns = self.introns[::-1]
            self.cdsintrons = self.cdsintrons[::-1]
        
if __name__ == "__main__":
    fin = readFile(args.input)
    fout = writeFile(args.output)

    geneTree = OrderedDict()
    for line in fin:
        if line.startswith("#"):
            continue
        tmp = line.strip().split("\t")
        if tmp[2] == "gene" or tmp[2] == "Selenocysteine" or tmp[2] == "start_codon" or tmp[2] == "stop_codon" or tmp[2] == "pseudogene": 
            continue
        tmp[3] = int(tmp[3])-1
        tmp[4] = int(tmp[4])
        d = parseStringgtf(tmp[8])
        if tmp[2] == "transcript":
            if d["gene_id"] not in geneTree:
                geneTree[d["gene_id"]] = gene(d["gene_id"],tmp[0],tmp[3],tmp[4],tmp[6])
            else:
                if tmp[0] != geneTree[d["gene_id"]].chrom or tmp[6] != geneTree[d["gene_id"]].strand:
                    sys.exit("{0} is error".format(d["gene_id"]))
            if int(tmp[3]) < geneTree[d["gene_id"]].start:
                geneTree[d["gene_id"]].start = int(tmp[3])
            if int(tmp[4]) > geneTree[d["gene_id"]].end:
                geneTree[d["gene_id"]].end = int(tmp[4])
            
            if d["transcript_id"] not in geneTree[d["gene_id"]].rna:
                geneTree[d["gene_id"]].rna[d["transcript_id"]] = trans(d["transcript_id"],tmp[0],tmp[3],tmp[4],tmp[6])
            else:
                sys.exit("there are two same transcript {0}".format(d["transcript_id"]))
        elif tmp[2] == "CDS":
            if d["transcript_id"] not in geneTree[d["gene_id"]].rna:
                sys.exit("you gtf file has wrong format1")
            geneTree[d["gene_id"]].rna[d["transcript_id"]].cds.append([tmp[3],tmp[4]])
        elif tmp[2] == "UTR":
            pass
        elif tmp[2] == "exon":
            if d["transcript_id"] not in geneTree[d["gene_id"]].rna:
                print(line)
                sys.exit("you gtf file has wrong format2")
            geneTree[d["gene_id"]].rna[d["transcript_id"]].exons.append([tmp[3],tmp[4]])
        else:
            sys.stderr.write(line)
            sys.exit("error2!\n")

    print("total gene: {0}".format(len(geneTree.keys())))
    for k in geneTree.keys():
        ge = geneTree[k]
        ge.check()
        fout.write("{0}\t{1}\t{2}\t{3}\tgene\t{4}\n".format(ge.chrom,ge.start,ge.end,ge.id,ge.strand))
        for j in ge.rna.keys():
            ts = ge.rna[j]
            ts.getIntron()
            ts.check()
            ts.getForword()
            fout.write("{0}\t{1}\t{2}\t{3}\ttranscript\t{4}\t{5}\n".format(ts.chrom,ts.start,ts.end,ge.id,ts.strand,ts.id))
            # for index,(sss,eee) in enumerate(ts.exons):
            #     fout.write("{0}\t{1}\t{2}\t{3}\texon\t{4}\t{5}\t{6}\n".format(ts.chrom,sss,eee,ge.id,ts.strand,ts.id,index))
            # if len(ts.introns) > 0:
            #     for index,[sss,eee] in enumerate(ts.introns):
            #         fout.write("{0}\t{1}\t{2}\t{3}\tintron\t{4}\t{5}\t{6}\n".format(ts.chrom,sss,eee,ge.id,ts.strand,ts.id,index))
            # if len(ts.cds) > 0:
            #     for index,(sss,eee) in enumerate(ts.cds):
            #         fout.write("{0}\t{1}\t{2}\t{3}\tcd\t{4}\t{5}\t{6}\n".format(ts.chrom,sss,eee,ge.id,ts.strand,ts.id,index))
            # if len(ts.cdsintrons) > 0:
            #     for index,[sss,eee] in enumerate(ts.cdsintrons):
            #         fout.write("{0}\t{1}\t{2}\t{3}\tcdsintron\t{4}\t{5}\t{6}\n".format(ts.chrom,sss,eee,ge.id,ts.strand,ts.id,index))

            ## ################################################################################### (刚转录出来的RNA) preRNA
            ## ######-------########---------#######------#########------------------------####### (成熟RNA) mRNA remove intron
            ##                               #######------#########                                (翻译区)

            if len(ts.exons) == 1:
                fout.write("{0}\t{1}\t{2}\t{3}\texon\t{4}\t{5}\t{6}\t{7}\n".format(ts.chrom,ts.exons[0][0],ts.exons[0][1],ge.id,ts.strand,ts.id,0,0))
                if len(ts.introns) > 0:
                    sys.exit("error! {0} has contradiction exon and intron".format(ts.id))
            elif len(ts.exons) > 1:
                if len(ts.introns) + 1 != len(ts.exons):
                    sys.exit("error! {0} has contradiction exon and intron".format(ts.id))
                exontss = 0
                tss = 0
                fout.write("{0}\t{1}\t{2}\t{3}\texon\t{4}\t{5}\t{6}\t{7}\n".format(ts.chrom,ts.exons[0][0],ts.exons[0][1],ge.id,ts.strand,ts.id,0,0))
                index = len(ts.introns)
                exontss += ts.exons[0][1]-ts.exons[0][0]
                tss += ts.exons[0][1]-ts.exons[0][0]
                for iiii in range(index):
                    fout.write("{0}\t{1}\t{2}\t{3}\tintron\t{4}\t{5}\t{6}\t{7}\n".format(ts.chrom,ts.introns[iiii][0],ts.introns[iiii][1],ge.id,ts.strand,ts.id,tss,iiii))
                    fout.write("{0}\t{1}\t{2}\t{3}\texon\t{4}\t{5}\t{6}\t{7}\n".format(ts.chrom,ts.exons[iiii+1][0],ts.exons[iiii+1][1],ge.id,ts.strand,ts.id,exontss,iiii+1))
                    exontss += ts.exons[iiii+1][1]-ts.exons[iiii+1][0]
                    tss += ts.introns[iiii][1]-ts.introns[iiii][0]
                    tss += ts.exons[iiii+1][1]-ts.exons[iiii+1][0]
            else:
                sys.exit("error! {0} has no exons!".format(ts.id))
    print("successful")
    fin.close()
    fout.close()
