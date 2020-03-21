#!/usr/bin/env python

# Program Header
# Course: Bi624
# Name:   Gerardo Perez
# Description: Our goal with this assignment
# is to create python script that visualizes motifs on sequences.
# The script also outputs the number of motifs for eache gene.
#
#
#
# Motif_mark
#
#
# Development Environment: Atom 1.38.2
# Version: Python 3.7.3
# Date:  03/14/2020
#################################################

# Programs Source Statements

# Imports module
import argparse
import re
import cairo
import math

# Creates an arguement passer
parser = argparse.ArgumentParser(description="A program that finds and colors motifs from genes. Max is 10 motifs")

# Adds arguemets by calling the arguement passer, for fasta file name
parser.add_argument("-f", "--fasta_file", help="Specify the fasta file with with genes", required=True)

# Adds arguemets by calling the arguement passer, for motif file name
parser.add_argument("-m", "--motifs", help="Specify the fasta file with with motifs", required=True)

# Parses arguemets through using the parse args method.
args = parser.parse_args()

# creates variables for args parse for input files
f=open(args.fasta_file,"r")
motif=open(args.motifs,"r")

# Stores the lines from files into variable.
lines=f.readlines()
motifs_lines=motif.readlines()

# colors for pycairo, 10 colors
col_list=['set_source_rgb(1, 0, 0)', #red
        'set_source_rgb(0, 0, 1)',   # royal blue
        'set_source_rgb(0.4, 0.9, 0.4)', #light green
        'set_source_rgb(0, 1, 1)',        #turquoise
        'set_source_rgb(1, 1, 0)',       # yellow
        'set_source_rgb(0.9, 0.4, 0.9)', # violet
        'set_source_rgb(0.2, 0.4, 0.8)', #college navy
        'set_source_rgb(0.3, 0.8, 0.6)', #cyan
        'set_source_rgb(1, 0.8, 0.5)',   #orange
        'set_source_rgb(1, 0, 1)'        # pink
        ]

# creates an empty dictionary
motifs_col={}

# valid bases found in motifs
valid_bases = ['y', 'g', 'c', 't', 'u']

# creates dictionary: key is valid base in motif, value is regex pattern to find base
comp_dict = {'y': '[c,t]', 'g': '[g]', 't': '[t]', 'c': '[c]', 'u':'[u,t]', 'a':'[a]'}

# Creates a function to get motifs from file to be stored as a key in a regex patter and a color as a value in a dictionary.
def get_motifs(motifs_lines):
    """Takes readlines in from motif file. returns a dictionary: keys are regex pattern and values are pycairo color scheme """
    count=0
    # for loop to itreate through each motifs lines
    for i in motifs_lines:

        # emptry string to concanate regex pattern
        str_motif=""

        #  removes new line character and splits at each tab., gets the 0th column (IMPORTANT!)
        i=i.strip().split("\t")[0]

        #for loop that iterates each charcter in line
        for j in i:

            # changes charecter to lower case
            j=j.lower()

            # adds the correspoding regex pattern of a charracter to a string
            str_motif=str_motif+comp_dict[j]

        # adds to dictionary: regex patter is key, pycairo color scheme are values
        motifs_col[str_motif]=col_list[count]

        count=count+1


    return(motifs_col)

# Creates a function to parse gene to gene name as key and the genes introns plus exons as values
def parse_gene(gene):
    """Takes readlines from fasta file. returns a dictionary: keys are the gene name header and values are its intron, exon, intron """

    # creates variable for function
    fasta=""
    dict={}
    mydict={}

    # for loop that will find the line with header and add the @ symbol at the end
    for i in gene:

        # check if character is in the line
        if ">" in i:

            # removes the new line character
            i=i.strip()

            # adds the symbol at the end of the line
            i=i+"@"

        # removes new line character and splits at each tab., gets the 0th column (IMPORTANT!)
        x=i.strip().split('\t')[0]

        # concanates a string for each line iterated
        fasta=fasta + x

    # splits the string at the carrot symbol into a list, store into a new variable
    genes=fasta.split(">")

    # removes the first index of the list
    genes=genes[1:len(genes)]

    # for loop that iterates through a string
    for i in genes:

        # splits the string at the @ symbol into a list store into a new
        i=i.split("@")

        # store keys as the gene name header and the values as the gene sequence.
        dict[i[0]]=i[1]

    # for loop that itereates through a dictionary keys as the gene name header and the values as the gene sequence.
    for k, v in dict.items():

        # splits gene sequence by first charcter uppercase
        split_seq=re.sub( r"([A-Z]+)", r" \1", v).split()

        mylist=[]

        # for loop that itereates the sequence that was split
        for i in split_seq:

            # checks if character is lower case
            if i.islower()==False:

                # if false split sequence  by first charcter lowercase
                sec_split=re.sub( r"([a-z]+)", r" \1", i).split()

                # for loop that itereates through splitted sequence
                for j in sec_split:

                    # adds sequence to list
                    mylist.append(j)

            # if character is lowerecase true, adds sequence to list
            else: mylist.append(i)

            # stores in a dictionary: keys are the gene name header and values are its intron, exon, intron
            mydict[k]=mylist

    return(mydict)

# Creates a function that finds gene locations then prints motif, gene header name and motif counts for each gene
def motifs_location(motifs, gene):
    """Takes readlines from fasta file and from motif file, prints motif, gene header name and motif counts for each gene"""

    motif_dict={}
    total_dict={}
    motif_list=[]
    count_motif=0


    # for loop that itereates through dictionary: keys are motif regex pattern and values are pycairo color scheme
    for key, v in get_motifs(motifs).items():

        # counter to  give space draw  the next gene
        counter=0
        color_scheme=v

        # compiles motif regex pattern for finditer() function
        p = re.compile(key)


        # for loop that itereates through dictionary: gene name as key and the genes introns plus exons as values
        for k, v in parse_gene(gene).items():

            # variable to be used in the function
            start_loc=[]
            end_loc=[]
            str=""
            loc_list=[]
            total=0

            # for loop that itereates through each gene intron, exon, intron
            for i in v:

                # concanates the intron/exon to a string
                str=str+i

            # for loop that finds motif regex pattern location,converts to lower case
            for m in p.finditer(str.lower()):

                # stores starting loction of motif to a list
                start_loc.append(m.start())

                # stores end loction of motif to a list
                end_loc.append(m.end())

            # counter to call the end location
            count=0

            # for loop that itereates through each motif starting location
            for i in start_loc:

                context.set_line_width(15)

                # sets a specific pycairo color scheme from color list
                get_col="context."+color_scheme

                # executes the string
                exec(get_col)

                # starts the line drawing at a specific location
                context.move_to(50+i, 150+counter)

                # stops the line drawing at a specific location
                context.line_to(50+end_loc[count],150+counter)
                context.stroke()

                count=count+1

            counter=counter+100

            # prints a specific motif
            print(motifs[count_motif].strip())

            # prints the gene name header and counts of the specific motif
            print("{0}: {1}".format(k,len(start_loc)))

        # increments by one to print the next motif
        count_motif=count_motif+1
        print("\n")
    return(motif_list)

# Creates a function that finds gene length, then returns dictionary with keys as gene headers and values as total gene length
def gene_length(motifs, gene):
    """Takes readlines from fasta file and from motif file, returns dictionary with keys as gene headers and values as total gene length"""

    # variable to be used in the function
    motif_dict={}
    total_dict={}

    # for loop that itereates through dictionary: keys are motif regex pattern and values are pycairo color scheme
    for k, v in get_motifs(motifs).items():


        color_scheme=v

        # compiles motif regex pattern for finditer() function
        p = re.compile(k)

        # for loop that itereates through dictionary: gene name as key and the genes introns plus exons as values
        for k, v in parse_gene(gene).items():

            # variable to be used in the function
            loc_list=[]
            total=0

            # for loop that itereates through each gene intron, exon, intron
            for i in v:

                # for loop that finds motif regex pattern location, converts to lower case
                for m in p.finditer(i.lower()):

                    # stores starting loction with the addition of the previous iterated gene length
                    loc_list.append(m.start()+total)

                # Sums the gene length after each gene iteration
                total=len(i)+total
            # stores in a dictionary: keys are the gene name header and values of total gene length
            total_dict[k]=total

    return(total_dict)





# Creates a function to draw the motif legend
def draw_motifs_legend(motifs):
    """Takes readlines from motif file to draw the motif alpha characters and its specific color"""
    space=0
    index=0
    counter_2=0

    # for loop that itereates through dictionary: keys are motif regex pattern and values are pycairo color scheme
    for k, v in get_motifs(motifs_lines).items():
        context.set_line_width(7)
        # Use the color scheme to print a specific color
        get_col="context."+v
        # coverts the string to a variable
        exec(get_col)
        context.move_to(50, 47+space)
        context.line_to(75,47+space)
        context.stroke()
        context.set_source_rgb(0, 0, 0)
        context.select_font_face("Arial bold", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(8)
        context.move_to(76, 49+space)
        context.show_text(motifs_lines[index])
        space=space+8
        index=1+index

# Creates a function to draw gene headers and gene lengths as lines
def draw_gene_len(motifs, gene):
    """Takes readlines from fasta file and from motif file to draw each gene header as text and the length of each gene as a line"""
    counter=0

    # for loop that itereates through dictionary: keys as gene headers and values as total gene length
    for k, v in gene_length(motifs, gene).items():
        context.set_source_rgb(0, 0, 0)
        context.select_font_face("Arial bold", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(10)
        context.move_to(50, 135+counter)
        # draws gene header as text
        context.show_text(k)
        context.set_line_width(1)
        context.set_source_rgb(0, 0, 0)
        context.move_to(50, 150+counter)
        # draws the length of the line
        context.line_to(50+v,150+counter)
        context.stroke()
        counter=counter+100

# Creates a function to draw the exon for each gene as a thick line
def draw_exon(gene):
    """Takes readlines from fasta file to draw the exon for each gene as a thick line"""
    counter=0
    total=0
    finish=0

    # for loop that itereates through dictionary: gene name as key and the genes introns plus exons as values
    for gene_name, int_ex in parse_gene(gene).items():
        start=0
        for i in int_ex:
            if i.islower()==False:
                #width bigger than gene length line
                context.set_line_width(15)
                context.set_source_rgb(0, 0, 0)
                context.move_to(50+start, 150+counter)
                context.line_to(50+start+len(i),150+counter)
                context.stroke()
                counter=counter+100
            start=start+len(i)



# parameters for device-space units
WIDTH, HEIGHT = 1000, 1000

# abstract type representing all different drawing targets
surface = cairo.SVGSurface("plot.svg", WIDTH, HEIGHT)

#actual drawings are performed
context = cairo.Context(surface)


#Calls the function to run the script
draw_exon(lines)
draw_gene_len(motifs_lines, lines)
draw_motifs_legend(motifs_lines)
motifs_location(motifs_lines, lines)


# Program Output (Commented out)
"""

ygcy
INSR chr19:7150261-7150808 (reverse complement): 12
ygcy
MBNL chr3:152446461-152447003: 10
ygcy
ATP2A1 chr16:28903467-28904044: 26
ygcy
CLASP1 chr2:121444593-121445363 (reverse complement): 18


GCAUG
INSR chr19:7150261-7150808 (reverse complement): 2
GCAUG
MBNL chr3:152446461-152447003: 2
GCAUG
ATP2A1 chr16:28903467-28904044: 0
GCAUG
CLASP1 chr2:121444593-121445363 (reverse complement): 1


catag
INSR chr19:7150261-7150808 (reverse complement): 1
catag
MBNL chr3:152446461-152447003: 1
catag
ATP2A1 chr16:28903467-28904044: 1
catag
CLASP1 chr2:121444593-121445363 (reverse complement): 0


YYYYYYYYYY
INSR chr19:7150261-7150808 (reverse complement): 1
YYYYYYYYYY
MBNL chr3:152446461-152447003: 2
YYYYYYYYYY
ATP2A1 chr16:28903467-28904044: 3
YYYYYYYYYY
CLASP1 chr2:121444593-121445363 (reverse complement): 6

"""
