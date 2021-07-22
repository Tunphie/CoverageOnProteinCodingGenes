#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import re
import pysam

description = '''
    Splits mapped reads in BAM format into their localization in exonic/intronic regions and UTRs.
   '''

epilog = "Author: Ann-Sophie Seistrup, as@seistrup.dk"

#######


def getArgs():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )


    parser.add_argument(
        '-gtf',
        required=True,
        type=str,
        help='gtf-file with annotation'
    )
    parser.add_argument(
        '-genes',
        required=True,
        type=str,
        help='genes to extract information for'
    )
    parser.add_argument(
        '-bam',
        nargs="+",
        required=True,
        type=str,
        help='bamfile with reads'
    )

    parser.add_argument(
        '-o','--out',
        required=False,
        type=str,
        default=".",
        help='Output directory'
    )

    args = parser.parse_args()
    return args

######


def readgtf(gtffile):
    """
        Reads a GTF file with the format:
        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
        Returns a dict with gene --> feature --> range (where feature is also "chr" )
    """

    gene_dict = {}

    relevant_feats = ["exon","gene","five_prime_utr","three_prime_utr"]

    with open(gtffile) as f:
        for line in f:
            fields = line.rstrip().split('\t')

            if len(fields) >= 9:  # looks like a proper GTF line and not a header line
                chr = str(fields[0])
                feature = str(fields[2])
                start = int(fields[3])
                end = int(fields[4])
                gene = re.sub('";.*', '', re.sub('.*gene_id "', '', str(fields[8])))

                if feature not in relevant_feats:
                #    print("no ", feature)
                    continue
                #else:
                #    print("yes",feature)

                if gene not in gene_dict.keys():
                    gene_dict[gene] = {}
                    gene_dict[gene]["chr"] = chr

                if feature not in gene_dict[gene].keys():
                    gene_dict[gene][feature] = [(start,end)]
                else:
                    gene_dict[gene][feature] += [(start, end)]
    #print(gene_dict)

    return gene_dict


def readgenes(genesfile):
    '''
        Returns alist of all gene IDs on separate lines of a .txt file
    '''
    genes = []
    with open(genesfile,"r") as f:
        for line in f:
            genes += [str(line.rstrip())]
    return genes


def traversebam(bamfile,genes,gtf):
    '''
        goes over all reads in a bamfile
    '''
    bam = pysam.AlignmentFile(bamfile)
    counts_dict = {}
    for gene in genes:
        if gene not in gtf.keys():
            print("Gene " + gene + " not found in gtf-file.")
            continue

        featlist = [0,0,0,0]
        featlistnames = ["exon","three_prime_utr", "five_prime_utr", "intron"]

        for read in bam.fetch(reference=gtf[gene]["chr"], start=gtf[gene]["gene"][0][0], end=gtf[gene]["gene"][0][1]):
            positions = read.get_reference_positions()
            start = positions[0]
            end = positions[-1]
            #print(start, end, gene)

            partial = 0
            found = False
            deduct = []
            for i in range(0,3):
                try:
                    found_in_rng = False
                    for rng in gtf[gene][featlistnames[i]]:
                        #print(rng, featlistnames[i])
                        if start >= rng[0] and end <= rng[1]: #Completely in feature
                            found = True
                            found_in_rng = True
                            #break
                            #featlist[i] +=1
                        elif start <= rng[0] and end >= rng[0] or start <= rng[1] and end >= rng[1]: # in feature but also in another feature
                            found = True
                            found_in_rng = True
                            partial += 1
                            deduct.append(i)
                            #featlist[i] += 1
                    if found_in_rng:
                        featlist[i] +=1
                except KeyError:
                    #print(gene, featlistnames[i])
                    continue
            if not found:
                featlist[3] += 1
            if partial == 1:
                if 0 in deduct:    #only count it to overlap with intron if it's exon. UTRs will be asusmed to be wrongly annotated. 
                    featlist[3] += 1
                    deduct.append(3)
                else:
                    partial = 0
                    deduct = []
            if len(deduct) > 2:
                deduct = list(set(deduct))
                if len(deduct) == 1:
                    partial == 0
            if partial > 0:
                for i in deduct:
                    featlist[i] -= .5
        #    if gene == "WBGene00043992":
        #        print(gene,deduct,featlist)
            del(deduct,partial)
        counts_dict[gene] = featlist
    return counts_dict


def writeoutput(counts_dict,outfile):
    with open(outfile,"w") as f:
        f.write("Gene\tExon\t3'-UTR\t5'-UTR\tIntron")
        for gene in counts_dict.keys():
            f.write("\n" + gene + "\t" + "\t".join(str(round(count,2)) for count in counts_dict[gene]))



######
if __name__ == '__main__':
    args = getArgs()
    gtffile = args.gtf
    genesfile = args.genes
    bamfiles = args.bam
    outdir = args.out


    print("reading gtf..............................................")
    gtf = readgtf(gtffile)
    print("reading genesfile........................................")
    genes = readgenes(genesfile)
    print("traversing bamfiles......................................")
    for bamfile in bamfiles:
        print(".............." + re.sub(".*/", "", bamfile) + "..............")
        counts_dict = traversebam(bamfile,genes,gtf)
        print("writing output file......................................")
        outfile = outdir + "/" + re.sub(".*/", "",re.sub(".bam", "", bamfile)) + ".countsperfeature.txt"
        writeoutput(counts_dict,outfile)
