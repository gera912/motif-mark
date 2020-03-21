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


def motifs_location(motifs, gene):
    """Takes readlines from fasta file and from motif file, prints motif and motif counts for each gene"""

    motif_dict={}
    total_dict={}
    motif_list=[]
    count_motif=0


    
    for key, v in get_motifs(motifs).items():

        counter=0
        color_scheme=v
        p = re.compile(key)


        for k, v in parse_gene(gene).items():
            start_loc=[]
            end_loc=[]

            str=""
            loc_list=[]
            total=0
            for i in v:
                str=str+i

            for m in p.finditer(str.lower()):
                start_loc.append(m.start())
                end_loc.append(m.end())

            #break
            count=0
            for i in start_loc:
                context.set_line_width(15)
                get_col="context."+color_scheme
                exec(get_col)
                context.move_to(50+i, 150+counter)

                context.line_to(50+end_loc[count],150+counter)
                context.stroke()
                count=count+1
            counter=counter+100

            print(motifs[count_motif].strip())
            print("{0}: {1}".format(k,len(start_loc)))
        count_motif=count_motif+1
        print("\n")
    return(motif_list)


def gene_length(motifs, gene):
    motif_dict={}
    total_dict={}
    for k, v in get_motifs(motifs).items():
        #print(k)
        color_scheme=v
        p = re.compile(k)
        for k, v in parse_gene(gene).items():
            loc_list=[]
            total=0
            for i in v:
                for m in p.finditer(i.lower()):
                    loc_list.append(m.start()+total)


                total=len(i)+total
            total_dict[k]=total
            motif_dict[k]=loc_list



    return(total_dict)



gene_length(motifs_lines, lines)

WIDTH, HEIGHT = 1000, 1000

surface = cairo.SVGSurface("plot.svg", WIDTH, HEIGHT)
context = cairo.Context(surface)


motif_dict={}
total_dict={}



def draw_motifs_legend(motifs):
    space=0
    index=0
    counter_2=0

    for k, v in get_motifs(motifs_lines).items():

        context.set_line_width(7)
        get_col="context."+v

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


def draw_gene_len(motifs, gene):
    counter=0
    for k, v in gene_length(motifs, gene).items():

        context.set_source_rgb(0, 0, 0)
        context.select_font_face("Arial bold", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(10)
        context.move_to(50, 135+counter)
        context.show_text(k)

        context.set_line_width(1)
        context.set_source_rgb(0, 0, 0)
        context.move_to(50, 150+counter)
        context.line_to(50+v,150+counter)
        context.stroke()
        counter=counter+100

def draw_exon(gene):
    counter=0
    total=0

    finish=0
    for gene_name, int_ex in parse_gene(gene).items():
        start=0
        for i in int_ex:
            if i.islower()==False:
                context.set_line_width(15)
                context.set_source_rgb(0, 0, 0)
                context.move_to(50+start, 150+counter)
                context.line_to(50+start+len(i),150+counter)
                context.stroke()
                counter=counter+100
            start=start+len(i)



draw_exon(lines)

draw_gene_len(motifs_lines, lines)


draw_motifs_legend(motifs_lines)


motifs_location(motifs_lines, lines)
