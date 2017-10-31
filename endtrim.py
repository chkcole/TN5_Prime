#!/usr/bin/env python2.7
"""
    Charles Cole chkcole@ucsc.edu 10/16/2014
    The purpose of this program is to parse FASTA or FASTQ files and print them
    in FASTA or FASTQ format with the possibility of trimming the sequences
    within the files based upon their quality. The input files should be in
    either FASTA or FASTQ format. The FASTA files should follow the NCBI FASTA
    format specifications which could be found at
    http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml as of the date this
    program was created. The Quality file which is paired with the FASTA
    file should be of the exact same format as the FASTA file except the
    sequence characters are replaced with integers representing the phred
    quality score of each position within the DNA/RNA/Protein sequence
    seperated with spaces. The FASTA and Quality files should be created
    so that the sequence occur in the same order within each file. for
    example if a sequence with an ID of
    "seq1" is the first in the FASTA file, then the quality scores for
    sequence with ID "seq1" should be the first in the quality file as well.

    The FASTQ files should be formated in the such that the following items
    appear in the exact order they are shown: An ID line where the first
    character is an "@" symbol. immediatly after that symbol any non-whitespace
    characters((" \t\n\r\f\v") in python set notation)
    will be interpreted as being part of the ID. After the first whitespace
    character other than "\n" any characters will be interpreted as being part
    of the comments. After "\n" all characters belonging to either the Nucleic
    acid set("ACGTNUKSYMWRBDHV") or the Amino Acid Set
    ("APBQCRDSETFUGVHWIYKZLXM*N-") will be interpreted as beloning to the
    DNA/RNA/Protein sequence, all other characers will be ignored. It will
    keep adding characters to the sequnce in this pattern until a "+"
    character is encountered. anything between this "+" symbol and the next
    "\n" will be ignored. Once a new line begins, it reads the next
    "n" bytes from the file where "n" is the number of bytes read in sequence.
    Basically, it assumes equal sized quality sequences and DNA/RNA/Protein
    sequences.

    The commands you use are as follows and were taken straight from Kevin 
    Karplus's website
    (http://users.soe.ucsc.edu/~karplus/bme100/f14/Parsing.html)
    detailing the specifications for this program:
--min_qual [integer]
    The lowest quality value that can appear in the output. If a base of
    lower quality appears in an input sequence, the sequence is truncated 
    to remove all bases at that position or later.
--in_33 [filename]
    Input file in fastq format, using Phred+33 coding
--in_64 [filename]
    Input file in fastq format, using Phred+64 coding
--out_33 [filename]
    Output file in fastq format, using Phred+33 coding
--out_64 [filename]
    Output file in fastq format, using Phred+64 coding
--in_fasta [filename]
    Input file in fasta format. If used, there must also be a
    --in_qual argument.
--in_qual [filename]
    Input file in fasta quality format. If used, there must also be a
    --in_fasta argument.
--out_fasta [filename]
    Output file in fasta format.
--out_qual [filename]
    Output file in fasta quality format. If used, there must also be a
    --out_fasta argument. It is possible to have a --out_fasta without
    --out_qual, in which case the quality information is quietly discarded.
    This is a fairly common use case, when converting data for programs that
    can't use quality information.    
"""
from __future__ import print_function, division
import sys
import argparse
import itertools
import string
import re
import read_FASTQ
def print_fastq(read,sink,offset,minimum_Quality):
    """
    This function prints a read into fastq format in either phred 33
    or phred 64 and trims the printed read such that the read and
    quality scores only extend to the point before which the first
    score is less than the minimum quality accepts a tuple consisting
    of (string,string,list,string) The first string represents the ID,
    the second string represents the sequence, the list of integers
    represents the quality scores and the last string is any comments
    that may appear.

    accepts a writable file object representing the file you want to
    print to. accepts a int indicating the offset for th output accepts
    an int representing the minimum quality score which all sequence
    characters should have.
    """
    ID = read[0] #this string holds the ID
    sequence = read[1] #this string holds the sequence
    quality = read[2] #this list holds the quality scores
    comment = read[3] #this string holds the comments, if any exist
    trimmed_sequence = "" # this string holds the sequence after trimming
    trimmed_quality = "" #this string holds the converted quality scored
                         #after trimming
    for character,score in  zip(sequence,quality):
        if score < minimum_Quality:
            break; 
        else:
            trimmed_sequence += character
            trimmed_quality += chr(score + offset)
    sink.write("@")
    sink.write(ID)
    if len(comment) > 0:
        sink.write(comment)
    sink.write("\n")
    sink.write(trimmed_sequence)
    sink.write("\n")
    sink.write("+\n")
    sink.write(trimmed_quality)
    sink.write("\n")

def print_fasta(read,fasta_sink,qual_sink,minimum_Quality,print_quality_flag):
    """
    This function prints a read into FASTA and Quality format and trims the
    printed read such that the read and quality scores only
    extend to the point before which the first score is less than the minimum
    quality accepts a tuple consisting of (string,string,list,string)
    The first string represents the ID, the second string represents the
    sequence, the list of integers represents the quality scores and the
    last string is any comments that may appear.

    accepts a writable file object representing the file you want to print to
    the FASTA data to.
    accepts a second writable file object representing the file you want to
    print the quality data to.
    accepts an int representing the minimum quality score which all sequence
    characters should have.
    accepts an int indicating if quality scores should be printed.
    1 >= yes, anything else = no
    """    
    ID = read[0] #this string holds the ID
    sequence = read[1] #this string holds the sequence
    quality = read[2] #this list holds the quality scores 
    comment = read[3] #this string holds the comments, if any exist
    trimmed_sequence = "" # this string holds the sequence after trimming
    trimmed_quality = [] # this list holds the quality scores after trimming
    for character,score in zip(sequence,quality):
        if score < minimum_Quality:
            break;
        else:
            trimmed_sequence += character
            trimmed_quality += [score]
    trimmed_quality = list(map(str,trimmed_quality))
    qual_string = " ".join(trimmed_quality)
    fasta_sink.write(">")
    fasta_sink.write(ID)
    if len(comment) > 0:
        fasta_sink.write(comment)
    fasta_sink.write("\n")
    fasta_sink.write(trimmed_sequence)
    fasta_sink.write("\n")
    if print_quality_flag >= 1: #only reports quality scored if there is a
                                #file to report them too
        qual_sink.write(">")
        qual_sink.write(ID)
        if len(comment) > 0:
            qual_sink.write(comment)
        qual_sink.write("\n")
        qual_sink.write(qual_string)
        qual_sink.write("\n")

def grab_Format_Parameters():
    input_file_list = [] # a list containing the file handles for the input
    output_file_list = [] # a list containing the file handles for the output
    input_offset = None # the offset to use when interpreting the fastq input 
    output_offset = None # the offset to use when encoding the fastq output
    minimum_quality = None##this integer indicates the minimum phred score a
                          ##character in the sequence
                          ##should have before trimming is activated
    
    param_Parser = argparse.ArgumentParser(description=__doc__)
    exclusive_phred_in = param_Parser.add_mutually_exclusive_group()
    exclusive_phred_out = param_Parser.add_mutually_exclusive_group()
    param_Parser.add_argument("--min_qual", help="Indicates the minimum \
    phred score before trimming is activated",required=True, type=int, \
    default= 0)
    exclusive_phred_in.add_argument("--in_33",help="Indicates that the \
    input fastq is in phred+33", \
    type=argparse.FileType("r"))
    exclusive_phred_in.add_argument("--in_64",help="Indicates that the \
    input fastq is in phred+64",type=argparse.FileType("r"))
    exclusive_phred_out.add_argument("--out_33",help="Indicates that the \
    output fastq is in phred+33",type=argparse.FileType("w"))
    exclusive_phred_out.add_argument("--out_64",help="Indicates that \
    the output fastq is in phred+64",type=argparse.FileType("w"))
    param_Parser.add_argument("--in_fasta",help="the file name of the input \
    fasta file. must be used with --in_qual [filename]", \
    type=argparse.FileType("r"))
    param_Parser.add_argument("--out_fasta",help="the file name of the  \
    output fasta file.",type=argparse.FileType("w"))
    param_Parser.add_argument("--in_qual",help="the file name of the \
    input quality file. must be used with --in_fasta [filename]", \
    type=argparse.FileType("r"))
    param_Parser.add_argument("--out_qual",help="the file name \
    of the output quality. must be used with --out_fasta [filename]", \
    type=argparse.FileType("w"))
    args = param_Parser.parse_args()

    if args.in_33 is not None:
        input_file_list.append(args.in_33)
        input_offset = 33

    if args.in_64 is not None:
        input_file_list.append(args.in_64)
        input_offset = 64

    if args.out_33 is not None:
        output_file_list.append(args.out_33)
        output_offset = 33

    if args.out_64 is not None:
        output_file_list.append(args.out_64)
        output_offset = 64

    if (args.in_fasta is not None) & (args.in_qual is not None) & \
    (input_offset is None):
        input_file_list.append(args.in_fasta)
        input_file_list.append(args.in_qual)
    elif (input_offset is None) | ((args.in_fasta is not None) & \
    (args.in_qual is not None) & (input_offset is not None)):
        sys.stderr.write("ERROR: input must be specified with either \
    --in_fasta [filename] AND --in_qual [filename] OR --in_fastq [filename]")
        param_Parser.print_usage(file=sys.stderr)
        exit()
        
    if (args.out_fasta is not None) & (output_offset is None):
        output_file_list.append(args.out_fasta)
    elif (output_offset is None) | ((args.out_fasta is not None) & \
    (input_offset is not None)):
        sys.stderr.write("ERROR: output must be specified with either \
    --out_fasta [filename] AND/OR --out_qual [filename] OR --out_fastq \
    [filename]")
        param_Parser.print_usage(file=sys.stderr)
        exit()
        
    if (args.out_qual is not None) & (args.out_fasta is not None) \
    & (output_offset is None):
        output_file_list.append(args.out_qual)
    elif ((args.out_fasta is not None) & (args.out_qual is not None) \
    & (input_offset is not None)) | ((args.out_qual is not None) & \
    (args.out_fasta is None)):
        sys.stderr.write("ERROR: output must be specified with either \
    --out_fasta [filename] AND/OR --out_qual [filename] OR --out_fastq \
    [filename]")
        param_Parser.print_usage(file=sys.stderr)
        exit()

    return(input_file_list,output_file_list,input_offset,output_offset, \
    args.min_qual)
                               
def main(args):
    input_file_list,output_file_list,input_offset,output_offset, \
    minimum_quality = grab_Format_Parameters()

    if input_offset is not None:
        for read in read_FASTQ.read_fastq_with_convert(input_file_list[0], input_offset):
            if output_offset:
                print_fastq(read,output_file_list[0],output_offset, \
                minimum_quality)
            else:
                if len(output_file_list) == 1:
                    print_fasta(read,output_file_list[0],"",minimum_quality,0)
                else:
                    print_fasta(read,output_file_list[0],output_file_list[1], \
                    minimum_quality,1)
                
    else:
        for read in read_FASTQ.read_fasta_with_quality(input_file_list[0], \
        input_file_list[1]):
            if output_offset:
                print_fastq(read,output_file_list[0],output_offset, \
                minimum_quality)
            else:
                if len(output_file_list) == 1:
                    print_fasta(read,output_file_list[0],"",minimum_quality,0)
                else:
                    print_fasta(read,output_file_list[0],output_file_list[1], \
                    minimum_quality,1)


if __name__ == "__main__" :
    sys.exit(main(sys.argv))
