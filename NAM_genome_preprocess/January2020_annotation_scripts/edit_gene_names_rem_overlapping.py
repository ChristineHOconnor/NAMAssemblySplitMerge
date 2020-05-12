"""This script will take a gff file and a file with current and new gene names and replace the old gene names in the gff file. A new gff file will be created and saved. Regular expressions are used to find and replace gene names. Originally I used python's replace function but that caused isssues becase searching for g7 would lead to replacements for g7 and also g70 etc. 
This version of the script will remove the 'outer' (longest) gene in any cases of nested genes
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
            -i or in_list list of genes that will be kept
            -l or reduced_list list of genes that will be excluded
        \n""")
#add in 2 arguements to get lists of genes that do and don't overlap with 
try:
    opts, args = getopt.getopt(sys.argv[1:], "g:o:r:k:i:l:h", ["gff_file=",
                                                         "output_file=",
                                                         "red_output_file=",
                                                         "name_key=",
                                                         "in_line=",
                                                         "reduced_list=",
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
    elif opt in ("-i", "--in_list"):
        in_genes = arg
    elif opt in ("-l", "--reduced_list"):
        red_list = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"

#replacing gene names works well and is quick
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
#02/13/2020: this is removing way too many genes 
#empty list to save the names of genes nested in other genes 
genes_for_final_gff = list()

prev_gene_end = 0
prev_gene_strand = ""
prev_gene_chr = "chr1"
prev_gene_type = ""
prev_gene_line = "first line"
#add count variable to control print statements
count = 1
n = 1
#find genes that are nested in other genes
with open(output_file, "r") as output,\
     open(in_genes, "w") as kept_genes,\
     open(red_list, "w") as exclude_genes:
    for line in output:
        if line.startswith("#"): 
            continue
        else:
            fields = line.split("\t")
            if fields[2] == "gene": 
                curr_gene_start = int(fields[3])
                curr_gene_strand = fields[6]
                curr_gene_chr = fields[0]
                #need to add check to avoid skipping over first line of 
                if curr_gene_start < prev_gene_end and curr_gene_strand == prev_gene_strand and curr_gene_chr == prev_gene_chr:
                    #if the start of the current gene is before the end of the previous gene, and they are on the same strand and chromosome they are overlapping
                    #overwrite the prev_gene information with the current gene information and move on
                    #exclude_genes.write(str(fields)) #this is wrong because the current gene is nested in the previous gene, it's the previous gene I want to remove
                    exclude_genes.write(str(prev_gene_line))
                    exclude_genes.write("\n")
                    #print(prev_gene_line)
                    #Now I want to overwrite the previous gene information because I don't want it to go anywhere
                    prev_gene_end = int(fields[4])
                    prev_gene_strand = fields[6]
                    prev_gene_chr = fields[0]
                    prev_gene_type = fields[1]
                    prev_gene_line = fields
                    if count < 15:
                        print("overlapping current gene start is: " + str(curr_gene_start))
                        print("previous gene end is: " + str(prev_gene_end))
                    count = count + 1
                else:
                    #if the start of the current gene is after the end of the previous gene, or strand and/or chromosome don't match
                    #1) save the name of the prev gene to the list of genes to save to final gff file
                    #if n < 15: 
                    #    print("non_overlapping gene: "+ fields[8])
                    #n = n + 1
                    kept_genes.write(str(prev_gene_line))
                    kept_genes.write("\n")
                    if prev_gene_type == "NAM": #getting some weird stuff here 
                        gene_name = prev_gene_line[8].split(":")[1].split(";")[0]
                        genes_for_final_gff.append(gene_name)
                    #    print("NAM gene: "+ gene_name)
                    elif prev_gene_type == "AUGUSTUS":
                        #try to not strip the new line character
                        gene_name = prev_gene_line[8].split("=")[1]
                        gene_name = gene_name.strip("\n")
                        #append the gene name to the list
                        genes_for_final_gff.append(gene_name)
                    #    print("AUGUSTUS gene: " + gene_name)
                    #2) update the prev_gene information with the current gene information
                    #So, don't save the gene name until it doesn't overlap with the next gene
                    prev_gene_end = int(fields[4])
                    prev_gene_strand = fields[6]
                    prev_gene_chr = fields[0]
                    prev_gene_type = fields[1]
                    prev_gene_line = fields #should this be fields or line? 
                    #print(prev_gene_line[8]) Assigning fields to prev_gene_line works as expected 
            else: 
                continue
#print(len(genes_for_final_gff))
#the length of this list is right for the number of genes that should be output
#it seems like something is going wrong here
#not all lines that should be output are
#will maybe need to do another regexp search 
#c = 1
#Now I have a list of gene should be written to the final gff file
#final loop to create new file with outer genes of nested genes excluded
with open(red_output_file, "w") as red_output,\
     open(output_file, "r") as output:
    #creating a dictionary won't work here becuase there's only one entry
    for line in output:
        if line.startswith("#") or 'assembly' in line:
            red_output.write(str(line))
            #red_output.write("\n")
            #is_gene_nested = "False" 
        else:
            Zm_gene_search = pattern_Zm.search(line)
            #if c < 10:
            #    print(Zm_gene_search)
            #    #print(Zm_gene_search.group(0))
            #c = c + 1
            abi_gene_search = pattern_abi.search(line)
            if Zm_gene_search and Zm_gene_search.group(0) in genes_for_final_gff:
                #if c < 24:
                #    print(Zm_gene_search)
                #    print(Zm_gene_search.group(0))
                red_output.write(line)
            elif abi_gene_search and abi_gene_search.group(0) in genes_for_final_gff:
                    #if c < 60: 
                    #    print(abi_gene_search)
                    #    print(abi_gene_search.group(0))
                    #    c = c + 1
                red_output.write(line)
            else: 
                continue

output.close()
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


#similar to what I want to do: 
#https://stackoverflow.com/questions/42982962/renaming-name-id-in-gffile

#end of script 
