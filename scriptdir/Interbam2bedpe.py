import argparse
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),".."))
from Util.Tools import readSam,writeFile

parser = argparse.ArgumentParser(description = " convert align file (sam or bam file) into bedpe file ")
parser.add_argument("-i","--input",required=True,help="the input sam/bam files.")
parser.add_argument("-o","--output",required=True,help="the output bedpe format file")
args = parser.parse_args()

if __name__ == "__main__":
    fout = writeFile(args.output)
    for infile in args.input.split(","):
        fin = readSam(infile)
        for read1 in fin:
            read2 = next(fin)
            info1 = read1.query_name
            info2 = read2.query_name
            if info1 != info2:
                sys.exit("reads name are not match")
            fout.write("{0}\t{1}\t{2}".format(read1.reference_name,read1.reference_start,read1.reference_end))
            strand = "+"
            if read1.is_reverse:
                strand = "-"
            fout.write("\t{0}".format(strand))
            if read1.get_tag("HG") != "":
                fout.write("\t{0}".format(read1.get_tag("HG")))
            else:
                fout.write("\t-")

            fout.write("\t{0}\t{1}\t{2}".format(read2.reference_name,read2.reference_start,read2.reference_end))
            strand = "+"
            if read2.is_reverse:
                strand = "-"
            fout.write("\t{0}".format(strand))
            if read2.get_tag("HG") != "":
                fout.write("\t{0}".format(read2.get_tag("HG")))
            else:
                fout.write("\t-")
            fout.write("\t{0}\n".format(info1))
        fin.close()
    fout.close()
