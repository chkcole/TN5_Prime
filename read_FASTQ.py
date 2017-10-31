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
#from scipy.sparse.linalg import eigsh
#from scipy.linalg import eigh
whitespace = set(" \t\n\r\f\v")#this set is used in the read_fasta
                               #function to determine where the
                               #ID ends in on the ID line
nucleic_string = "ACGTNUKSYMWRBDHV"
amino_string = "APBQCRDSETFUGVHWIYKZLXM*N-"
FASTA_Nucleic_Sequence_Character_Set = set(nucleic_string.lower()).union( \
set(nucleic_string))
FASTA_Amino_Acid_Sequence_Character_Set = set(amino_string.lower()).union( \
set(amino_string))

#this set determines which characters the read_fasta and
#read_fastq_generic functions will allow in the seqeuence string. 
combined_Amino_Nucleic_Character_Set = \
FASTA_Nucleic_Sequence_Character_Set.union( 
FASTA_Amino_Acid_Sequence_Character_Set)

#this set is used in the read_quality function to determine which characters
#it will allow in the quality string.
quality_Sequence_Character_Set = set("1234567890 \t.")

def read_Bytes_From_Source(source):
    """
    this generator function yields one byte from the source at a time
    until the end of the file.
    """    
    byte = source.read(1)#fetches a single character byte from the source

    while(byte):
        yield(byte)
        byte = source.read(1)
    return

def read_fasta_Generic(source, Sequence_Set,start_ID_Line = '>'):
    """
    The purpose of this function is to parse sequences from a file with
    a generalized FASTA format. This function accepts a file handle
    which contains the FASTA data, a set object which contains the legal
    characters for the sequences and a character which indicates
    the beginning of a sequence header. This function generates a tuple
    object of the form (string,string,string) where the first string
    represents the ID, the second string represents the sequence and
    the third string is any comments that may appear NOTE: the comment
    string may be empty source should be a readable file in either
    FASTA or Quality format Sequence_Set should be the set of legal
    characters which you would expect to be in that sequence, all other
    will be skipped. start_ID_Line is the character which you expect
    to deliniate the start of an ID line.
    """
    is_Comment_Section = 0 #the character is part of a comment
    is_Sequence_Section = 0 #the character is part of a sequence
    is_ID_Section = 0 #the character is part of an ID
    is_Skipped_section = 0 #the character is part of a section which should be
                           #skipped. this is invoked when a ';' character is
                           #encountered and ends when a new line starts
    is_not_empty = 0 #used to make sure something has been encountered in
                     #the file
    ID = "" #this string holds the ID
    comment = "" #this string holds the comment, if there is any.
    sequence = "" #this string holds the sequence
    for byte in read_Bytes_From_Source(source):
        if is_ID_Section:
            if byte not in whitespace: #whitespace cannot be included in an ID
                ID += byte
            elif byte != "\n": #if these is a whitespace character other
                               #than newline, start adding to comment
                is_Comment_Section = 1
                is_ID_Section = 0
            else: #if there's a newline character, start adding characters
                  #to the sequence
                is_Sequence_Section = 1
                is_ID_Section = 0
                
        elif is_Comment_Section:
            if byte != "\n": #any character will be added to the comment,
                             #as long as it's on one line
                comment += byte
            else: #once comment is done, start adding characters
                  #to sequence
                is_Comment_Section = 0
                is_Sequence_Section = 1
                
        elif is_Sequence_Section: 
            if byte in Sequence_Set: #if character is part of legal
                                     #set, add to sequence string
                sequence += byte
            elif byte == ";": # this character means there is comment
                              #within the sequence string
                is_Sequence_Section = 0
                is_Skipped_Section = 1
            elif byte == start_ID_Line: # looks for begining of header,
                                        #start new sequence
                yield(ID,sequence,comment.strip())
                ID = ""
                comment = ""
                sequence = ""
                is_ID_Section = 1
                is_Sequence_Section = 0
        elif is_Skipped_section: 
            if byte == "\n": # will ignore everything after ";" until
                             #new line begins
                is_Skipped_Section = 0    
                is_Sequence_Section = 1
        else:
            if byte == start_ID_Line: #initiates collection of sequences.
                                      #if no state exists
                is_ID_Section = 1
                is_not_empty = 1
            else:
                sys.stderr.write("ERROR: not FASTA format, " \
                +start_ID_Line+" was not first character \
                encountered\n")
                exit()
                
    if is_not_empty:
        yield(ID,sequence," "+comment)

def read_fasta(source):
    """
    This function generates read tuples from a file in FASTA format
    This function generates a tuple object of the form
    (string,string,string) where the first string represents the
    ID, the second string represents the sequence and the third
    string is any comments that may appear
    NOTE: the comment string may be empty
    """
    for sequence_Tuple in read_fasta_Generic(source, \
    combined_Amino_Nucleic_Character_Set):
        yield(sequence_Tuple)
    source.close()

    
def read_quality(source):
    """
    This function generates read tuples from a file in Quality format
    This function generates a tuple object of the form (string,list)
    where the first string represents the ID, and the list of integers
    represents the quality scores.
    """
    for sequence_Tuple in read_fasta_Generic(source, \
    quality_Sequence_Character_Set):
        qualities = sequence_Tuple[1].split()
        qualities = map(int,qualities)
        qualities = list(qualities)
        sequence_Tuple = (sequence_Tuple[0],qualities)
        yield(sequence_Tuple)
    source.close()

def read_fasta_with_quality(fasta_file,quality_file):
    """
    this function accepts two file names in FASTA and Qual format
    respectively The function returns a tupel consisting of
    (string,string,list,string) The first string represents the ID,
    the second string represents the sequence, the list of integers
    represents the quality scores and the last string is any comments
    that may appear.
    """
    fasta_iterator = iter(read_fasta(fasta_file))
    quality_iterator = iter(read_quality(quality_file))

    #checks both generators to see if there are still lines left to be
    #generated. I believe if one file has more sequences than the other,
    #it will only return a number of sequences equal to the smallest
    #number of sequences in either file
    while [fasta_iterator,quality_iterator]: 
        sequence = next(fasta_iterator)
        quality = next(quality_iterator)
        yield(sequence[0],sequence[1],quality[1],sequence[2])


def read_fastq_generic(source, Sequence_Set, start_ID_Line= '@', \
start_Quality_Line ="+", offset = 33):
    """
    The purpose of this function is to parse sequences from a file
    with a generalized FASTQ format. This function accepts a file
    object, a set object and a string object This function generates a tuple
    object of the form (string,string,string,string) where the first
    string represents the ID, the second string represents the sequence and
    the third string is the quality scores and the fourth string is any
    comments that may appear.
    NOTE: the comment string may be empty.
    source should be a readable file in either FASTQ format Sequence_Set
    should be the set of legal characters which you would expect to be in
    that sequence, all other will be skipped. start_ID_Line is the
    character which you expect to deliniate the start of an ID line.
    start_Quality_Line is the character which you would expect to
    deliniate the start of a quality line NOTE: This parser relies
    heavily on the fact the stirng of quality scores and sequence
    strings in the file are the exact same length if they are not,
    function will print an error and exit.
    """
    State_variable = 0 # holds state of machine
                          # 0 = beginning of file
                          # 1 = ID section
                          # 2 = sequence section
                          # 3 = quality section
                          # 4 = comment_section
                          # 5 = quality ID section
    ID = None #this string holds the ID
    comment = ""#this string holds the comment, if there is any.
    sequence = ""#this string holds the sequence
    quality = ""#this string holds the quality scores
    ID_Line = "" #holds all bytes processed from ID Line
    for byte in read_Bytes_From_Source(source):
        if State_variable == 0:
            if byte == start_ID_Line:
                State_variable = 1
            elif byte not in whitespace:
                sys.stderr.write("ERROR: not FASTQ format, " \
                +start_ID_Line+" was not first character encountered\n"+byte+"\n")
                exit()
        elif State_variable == 1:
            if byte != "\n":
                ID_Line += byte
            else:
                ID,comment = process_ID_Line(ID_Line)
                if ID is None:
                    ID = ""
                if comment is None:
                    comment = ""
                State_variable = 2
        elif State_variable == 2:
            if byte in Sequence_Set:#if character is part of legal
                                    #set, add to sequence string
                sequence += byte
            elif byte == ";":#this character means there is comment
                             #within the sequence string
                State_variable = 4
            elif byte == start_Quality_Line:
                State_variable = 5
        elif State_variable == 4:
            if byte == "\n":
                State_variable = 2
                
        elif State_variable == 5:
            if byte == "\n":
                State_variable = 3
                
        elif State_variable == 3:
            if len(quality) < len(sequence):
                if ((ord(byte) - offset) >= 0):
                    quality += byte
                else:
                    sys.stderr.write("WARNING: score in sequence"+ID+"was less \
		    than zero\n")
                    sys.stderr.write(str(ord(byte) - offset)+"\n")
                    sys.stderr.write(str(len(quality))+"\n")
            else:
                yield(ID,sequence,quality,comment)
                ID = None
                comment = ""
                sequence = ""
                quality = ""
                ID_Line = ""
                State_variable = 0
    if ID is not None:
        yield(ID,sequence,quality,comment)

def process_ID_Line(ID_Line):
    ID_Line = ID_Line.strip()
    whitespace_regex = re.compile("(\S+)?(\s.+)?")
    ID_List = whitespace_regex.match(ID_Line).groups()
    return ID_List


def read_fastq(source,offset = 33):
    """
    This function reads a file in FASTQ format and returns a
    tuple consisting of (string,string,string,string)
    the first string represents the ID, the second string
    represents the sequence and the third string is the quality scores
    and the fourth string is any comments that may appear.
    the file should be readable
    """
    for reads in read_fastq_generic(source, \
    combined_Amino_Nucleic_Character_Set,"@","+",offset):
        yield reads
    source.close()


def read_fastq_with_convert(file, offset = 33):
    """
    This function reads a file in FASTQ format and returns a tuple
    consisting of (string,string,list,string) The first string
    represents the ID, the second string represents the sequence,
    the list of integers represents the quality scores and the last
    string is any comments that may appear. the offset is an integer
    which will indicate the encoding of the quality string
    """
    for read in read_fastq(file,offset):
        scores = list(map(lambda x: ord(x) - offset ,read[2]))
        yield (read[0],read[1],scores,read[3])

