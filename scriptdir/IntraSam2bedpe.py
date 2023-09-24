import argparse
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import readSam,writeFile,getStrand

def checkPairtagRead(read1,read2):
    tmp1 = read1.query_name
    tmp2 = read2.query_name
    if tmp1 != tmp2:
        return True
    return False

parser = argparse.ArgumentParser(description = " convert align file (sam or bam file) into bedpe file ")
parser.add_argument("-i","--input",required=True,help="the input sam/bam files. comma separated (file1,file2,file3,...  )")
parser.add_argument("-o","--output",required=True,help="the output bedpe format file")
args = parser.parse_args()

if __name__ == "__main__":
    fout = writeFile(args.output)
    for infile in args.input.split(","):
        fin = readSam(infile)
        for read1 in fin:
            read2 = next(fin)
            if checkPairtagRead(read1,read2):
                sys.stderr.write("two reads are not match \n{0}\n{1}\n".format(read1.tostring(),read2.tostring()))
                sys.exit()
            info1 = read1.query_name
            info2 = read2.query_name
            if read1.is_reverse:
                strand1 = "-"
            else:
                strand1 = "+"
            if read2.is_reverse:
                strand2 = "-"
            else:
                strand2 = "+"
            if strand1 != strand2:
                sys.stderr.write("two reads are not match \n{0}\n{1}\n".format(read1.tostring(),read2.tostring()))
                sys.exit()
            HGset1 = set(read1.get_tag("HG").split(";"))
            HGset2 = set(read2.get_tag("HG").split(";"))
            HGset = HGset1.intersection(HGset2)
            if len(HGset) == 0:
                sys.stderr.write("two reads are not match \n{0}\n{1}\n".format(read1.tostring(),read2.tostring()))
                sys.exit()
            for kk in HGset:
                if read1.reference_start <= read2.reference_start:
                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                        read1.reference_name,read1.reference_start,read1.reference_end,read2.reference_name,read2.reference_start,read2.reference_end,info1,kk))
                else:
                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                        read2.reference_name,read2.reference_start,read2.reference_end,read1.reference_name,read1.reference_start,read1.reference_end,info1,kk))
        fin.close()
    fout.close()