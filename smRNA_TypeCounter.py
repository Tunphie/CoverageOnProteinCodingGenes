#!/usr/bin/env python
import argparse
import pysam

description = """
     Traverses BAM formatted files and counts different C. elegans smRNA types using strict requirements:
     22Gs must start with a G and be between 19 and 23 nt
     26Gs must be 26nt and start with a G
     piRNAs/21Us must be between 19 and 22 nt, and map to a provided reference. No 5' nucleotide is required, but sequences of this length starting with a G will be assumed to be a 22G with no sequence being allowed to be counted as twice as two different types. 
     miRNAs must map to a provided miRNA reference, be less than 27nt in length, and have not been defined by any of the previous classes.
     All RNAs longer than 27nt will also be counted. 
"""

parser = argparse.ArgumentParser()

parser.add_argument("bamfiles", nargs="+")

parser.add_argument('--mirna',
                    required=True,
                    help='Reference miRNA GTF file')
parser.add_argument('--pirna',
                    required=True,
                    help='Reference piRNA GTF file')
args = parser.parse_args()

def read_gtf(gtffile):
    with open(gtffile) as f:
        dict = {}
        for line in f:
            fields = line.rstrip().split('\t')
            if len(fields) >= 9: # looks like a proper GTF line and not a header line
                chr = str(fields[0])
                feature = str(fields[2])
                begin = int(fields[3])
                end = int(fields[4])
                if not chr in dict.keys():
                    dict[chr] = {}
                dict[chr][begin,end] = 1
    return dict



def readBam(bamfile,mirna_dict,pirna_dict):
    file = pysam.AlignmentFile(bamfile, "rb")
    dict = {"21U":0,"22G":0,">27":0,"26G":0,"miRNA":0, "all":0}
    for chr in ["I", "II", "III","IV","V","X", "MtDNA"]:
        for read in file.fetch(chr):
            if read.is_unmapped:
                next
            dict["all"] += 1
            len = read.query_alignment_length
            fiveprime = read.query_alignment_sequence[0]
            written = False
            if 19 < len < 24:
                if read.is_reverse and fiveprime == "G":
                    dict["22G"] += 1
                    written = True
            elif len == 26 and fiveprime == "G":
                dict["26G"] += 1
                written = True
            elif len > 27:
                dict[">27"] += 1
                written = True
            if not written and not read.is_reverse: 
                start = read.get_reference_positions()[0]
                end = start+len
                if 19 < len < 23:
                    try:
                        for seg in pirna_dict[chr].keys():
                            if start >= seg[0] and end <= seg[1]:
                                dict["21U"] += 1
                                written = True
                                break
                    except KeyError:
                        continue
                if not written: 
                        try:
                            for seg in mirna_dict[chr].keys():
                                if start >= seg[0] and end <= seg[1]:
                                    dict["miRNA"] += 1
                                    written = True
                                    break
                        except KeyError:
                            continue
    for key in dict.keys():
        dict[key] = str(dict[key])
    return dict

def write(dict,file):
    with open(file,"a") as f:
        for key in dict.keys():
            f.write("\t".join([key,dict[key]["21U"],dict[key]["22G"],dict[key]["26G"],dict[key][">27"],dict[key]["miRNA"],dict[key]["all"]]) )

def main():
    pirna_dict = read_gtf(args.pirna)
    mirna_dict = read_gtf(args.mirna)
    dict = {}
    for file in args.bamfiles:
        dict[file] = {}
        dict[file] = readBam(file,mirna_dict,pirna_dict)
    write(dict,"counts.txt")

if __name__ == '__main__':
    main()


