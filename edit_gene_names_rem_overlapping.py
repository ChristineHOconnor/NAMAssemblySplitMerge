"This script will take a gff file and a file with current and new gene names and replace the old gene names in the gff file. A new gff file will be created and saved. Regular expressions are used to find and replace gene names. Originally I used python's replace function but that caused isssues becase searching for g7 would lead to replacements for g7 and also g70 etc." 
"""This version of the script will save lists of genes that do and don't overlap with previous genes. This will allow me to save the size distribuiton of each category of gene and origin annotation. 
Output files: new gff3 with all genes, new gff3 with next genes removed, list of non-nested genes, list of nested genes (that will be removed)
"""

import sys
import getopt
import operator
import re
import tempfile

#usage etc setup taken from Meesh's script
def usage():
    print("""\n
        python3 edit_gene_name.py
            -g or gff_file gff file with original version of names
            -o or output_file  Name of new gff file that will be created
            -r or red_output_file Name of new gff file with the nested genes removed
            -k or name_key File with current and new gene names
            -n or non_overlap File with genes that don't overlap with other genes
            -f or file_overlap File with genes that do overlap with other genes
        \n""")
#add in 2 arguements to get lists of genes that do and don't overlap with 
try:
    opts, args = getopt.getopt(sys.argv[1:], "g:o:r:k:n:f:h", ["gff_file=",
                                                         "output_file=",
                                                         "red_output_file=",
                                                         "name_key=",
                                                         "non_overlap=",
                                                         "file_overlap=",
                                                         "help"])

except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-g", "--gff_file"):
        gff_file = arg
    elif opt in ("-o", "--output_file"):
        output_file = arg
    elif opt in ("-r", "--red_output_file"):
        red_output_file = arg
    elif opt in ("-k", "--name_key"):
        name_file = arg
    elif opt in ("-n", "--non_overlap"):
        not_overlap = arg
    elif opt in ("-f", "--file_overlap"):
        file_w_overlap = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"

#use regexp instead of replace
#pattern = re.compile(r'(?<=gene:)(.*?)(?=;)|(?<=Name=)(.*?)(?=_)')
#pattern_all = re.compile(r'(?<=gene:)(.*?)(?=;)|(?<=Name=)(.*?)(?=_)|(?<=transcript:)(.*?)(?=_)|(?<=CDS:)(.*?)(?=_)|(?<=transcript_id=)(.*?)(?=_)|(?<=exon_id=)(.*?)(?=_)|(?<=protein_id=)(.*?)(?=_)|(?<=transcript:)(.*?)(?=$)|(?<=ID=)(.*?)(?=\.)|(?<=ID=)(.*?)(?=$)|(?<=Parent=)(.*?)(?=\.)|(?<=Parent=)(.*?)(?=$)')
pattern_Zm = re.compile(r'(?<=gene:)(.*?)(?=;)|(?<=Name=)(.*?)(?=_)|(?<=transcript:)(.*?)(?=_)|(?<=CDS:)(.*?)(?=_)|(?<=transcript_id=)(.*?)(?=_)|(?<=exon_id=)(.*?)(?=_)|(?<=protein_id=)(.*?)(?=_)|(?<=transcript:)(.*?)(?=$)')
pattern_abi = re.compile(r'(?<=ID=)(.*?)(?=\.)|(?<=ID=)(.*?)(?=$)|(?<=Parent=)(.*?)(?=\.)|(?<=Parent=)(.*?)(?=$)')
#there's a problem with my pattern searches, so I'll add 1 at a time until I get an issue
#problem: I didn't escape '.' 
#new problem: Zm pattern's aren't being found 
#01/21/2020: I decided to create 2 pattern search's, 1 for Zm gene names and 1 for ab initio g* gene names
#this works on test data and on the full gff file
#I'm not sure why, but when I have all the regular expressions in 1 re.compile statement not everything is found 

#write new gene names to output file with all genes 
with open(name_file) as gene_names,\
     open(gff_file) as f,\
     open(output_file, "w") as output:
     # A dictionary mapping gene names to what they will be changed to
    rename_dict = dict(line.split() for line in gene_names)
    #print(rename_dict) #this is fine
    for line in f:
        if line.startswith("#") or 'assembly' in line:
            output.write(line)
        else: 
            # Search for the gene name
            result1 = pattern_Zm.search(line)
            result2 = pattern_abi.search(line)
            # If we found a gene name and we can rename it then we substitute it
            if result1 and result1.group(0) in rename_dict:
                #print("present")
                #replace instances of 'gene:' and 'transcript:' with nothing
                new_line = pattern_Zm.sub(rename_dict[result1.group(0)], line)
                output.write(new_line)
                #out_f.write(new_line)
            else: 
                #continue
                if result2 and result2.group(0) in rename_dict:
                    new_line = pattern_abi.sub(rename_dict[result2.group(0)], line)
                    output.write(new_line)
                    #out_f.write(new_line)
                else: 
                    continue
   
#all the gene names have been replaced 
#now I need to get the list of genes that are nested in other genes 
#empty list to save the names of genes nested in other genes 
nested_gene_list = list()

prev_gene_end = 0
prev_gene_strand = ""
prev_gene_chr = ""

#find genes that are nested in other genes
with open(not_overlap, "w") as non_o,\
     open(file_w_overlap, "w") as file_o,\
     open(output_file, "r") as output:
    for line in output:
        if line.startswith("#"): 
            continue
        else:
            fields = line.split("\t")
            if fields[2] == "gene": 
                curr_gene_start = int(fields[3])
                curr_gene_strand = fields[6]
                curr_gene_chr = fields[0]
                if curr_gene_start < prev_gene_end and curr_gene_strand == prev_gene_strand and curr_gene_chr == prev_gene_chr:
                    #if the start of the current gene is before the end of the previous gene, and they are on the same strand and chromosome they are overlapping
                    #print(curr_gene_start, prev_gene_end, fields[8], end = "")
                    file_o.write(line)
                    if fields[1] == "NAM": 
                        gene_name = fields[8].split(":")[1].split(";")[0]
                        nested_gene_list.append(gene_name)
                    elif fields[1] == "AUGUSTUS":
                        gene_name = fields[8].split("=")[1]
                        gene_name = gene_name.strip("\n")
                        #append the gene name to the list
                        nested_gene_list.append(gene_name)
                else: 
                    prev_gene_end = int(fields[4])
                    prev_gene_strand = fields[6]
                    prev_gene_chr = fields[0]
                    non_o.write(line)
            else: 
                continue
#Now I have a list of gene that are nested in other genes
#final loop to create new file with nested genes excluded
with open(red_output_file, "w") as red_output,\
     open(output_file, "r") as output:
    for line in output:
        is_gene_nested = "False" 
        for gene in nested_gene_list: 
            if gene in line: 
                is_gene_nested = "True"
            else: 
                continue
        if is_gene_nested == "False": 
            #if the gene name in the line never appears 
            red_output.write(line)

output.close()
file_o.close()
non_o.close()
red_output.close()

#now replace all instances of 'gene:' and 'transcript:' with nothing
#Why? Because extra characters between ID= or Parent= and start of gene name will interfere with other steps in the pipeline
with open(output_file, "r") as sources:
    lines = sources.readlines()
with open(output_file, "w") as sources:
    for line in lines:
        sources.write(re.sub(r'=gene:|=transcript:', '=', line))
        #sources.write(re.sub(r'transcript:', '', line)
sources.close()

with open(red_output_file, "r") as sources_2:
    entries = sources_2.readlines()
with open(red_output_file, "w") as sources_2: 
    for entry in entries: 
        sources_2.write(re.sub(r'=gene:|=transcript:', '=', entry))
sources_2.close()

#Why: genes that overlap with other genes cause issues witht he parse longest ts py script
#Also, genes nested in other genes are unlikely to be functional
#1)get list of genes that are nested in other genes
##open(output_file, "r") as output,\                                                                                                                                                                                                         
#     #open(not_overlap, "w") as non_o,\  
"""with open(not_overlap, "w") as non_o,\
     open(file_w_overlap, "w") as file_o:
    for line in temp:
        if line.startswith("#"): 
            continue
        else:
            fields = line.split("\t")
            if fields[2] == "gene": 
                curr_gene_start = int(fields[3])
                curr_gene_strand = fields[6]
                curr_gene_chr = fields[0]
                if curr_gene_start < prev_gene_end and curr_gene_strand == prev_gene_strand and curr_gene_chr == prev_gene_chr:
                    #print(curr_gene_start, prev_gene_end, fields[8], end = "")
                    file_o.write(line)
                    if fields[1] == "NAM": 
                        gene_name = fields[8].split(":")[1].split(";")[0]
                        nested_gene_list.append(gene_name)
                    elif fileds[1] == "AUGUSTUS":
                        gene_name = fields[8].split("=")[1]
                        gene_name = gene_name.strip("\n")
                        #append the gene name to the list
                        nested_gene_list.append(gene_name)
                else: 
                    prev_gene_end = int(fields[4])
                    prev_gene_strand = fields[6]
                    prev_gene_chr = fields[0]
                    non_o.write(line)
            else: 
                continue

file_o.close()
non_o.close()
#2)if the gene name is found in a line, remove that line
#loop through the temp file
#if the genes in the nested_gene_list never appear in the line, write the line to the output file
with open(output_file, "w") as output:
    for line in temp:
        is_gene_nested = "False" 
        for gene in nested_gene_list: 
            if gene in line: 
                is_gene_nested = "True"
            else: 
                continue
        if is_gene_nested == "False": 
            output.write(line)

#close the output file
#out_f.close()

#now replace all instances of 'gene:' and 'transcript:' with nothing
#Why? Because extra characters between ID= or Parent= and start of gene name will interfere with other steps in the pipeline
with open(output_file, "r") as sources: 
    lines = sources.readlines()
with open(output_file, "w") as sources:
    for line in lines: 
        sources.write(re.sub(r'=gene:|=transcript:', '=', line))
        #sources.write(re.sub(r'transcript:', '', line)
sources.close()
"""
#similar to what I want to do: 
#https://stackoverflow.com/questions/42982962/renaming-name-id-in-gffile

#end of script 
