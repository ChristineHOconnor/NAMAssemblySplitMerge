"This script is used to replace cds fasta gene names, and the canonical gene list names, with the new gene names. This script is basically the same thing as the edit gff script but the regular expressions are slightly different"

import sys
import getopt
import operator
import re

#usage etc setup taken from Meesh's script
def usage():
    print("""\n
        python3 edit_gene_name.py
            -f or fasta_file fasta file with original version of names
            -o or output_file  Name of new gff file that will be created 
            -k or name_key File with current and new gene names
        \n""")
try:
    opts, args = getopt.getopt(sys.argv[1:], "f:o:k:h", ["fasta_file=",
                                                         "output_file=",
                                                         "name_key="
                                                         "help"])

except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-f", "--fasta_file"):
        fasta_file = arg
    elif opt in ("-o", "--output_file"):
        output_file = arg
    elif opt in ("-k", "--name_key"):
        name_file = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"

#use regexp instead of replace
#pattern = re.compile(r'(?<=gene:)(.*?)(?=;)|(?<=Name=)(.*?)(?=_)')
#pattern_all = re.compile(r'(?<=gene:)(.*?)(?=;)|(?<=Name=)(.*?)(?=_)|(?<=transcript:)(.*?)(?=_)|(?<=CDS:)(.*?)(?=_)|(?<=transcript_id=)(.*?)(?=_)|(?<=exon_id=)(.*?)(?=_)|(?<=protein_id=)(.*?)(?=_)|(?<=transcript:)(.*?)(?=$)|(?<=ID=)(.*?)(?=\.)|(?<=ID=)(.*?)(?=$)|(?<=Parent=)(.*?)(?=\.)|(?<=Parent=)(.*?)(?=$)')
pattern_Zm = re.compile(r'(?<=\>)(.*?)(?=_)|(?<=^Z)(.*?)(?=_)')
pattern_abi = re.compile(r'(?<=\>)(.*?)(?=\.)|(?<=^g)(.*?)(?=\.)')
#2 pattern searches
#this works on test data and on the full gff file
#I'm not sure why, but when I have all the regular expressions in 1 re.compile statement not everything is found 

with open(name_file) as gene_names,\
     open(fasta_file) as f,\
     open(output_file, 'w') as out_f:
    # A dictionary mapping gene names to what they will be changed to
    rename_dict = dict(line.split() for line in gene_names)
    #print(rename_dict) #this is fine
    for line in f:
        if line.startswith(">"):
            result1 = pattern_Zm.search(line)
            result2 = pattern_abi.search(line)
            # If we found a gene name and we can rename it then we substitute it
            if result1 and result1.group(0) in rename_dict:
                #print("present")
                #replace instances of 'gene:' and 'transcript:' with nothing
                new_line = pattern_Zm.sub(rename_dict[result1.group(0)], line)
                out_f.write(new_line)
            else:
                #continue
                if result2 and result2.group(0) in rename_dict:
                    new_line = pattern_abi.sub(rename_dict[result2.group(0)], line)
                    out_f.write(new_line)
                else:
                    continue
        else: 
            out_f.write(line)

out_f.close()
f.close()
gene_names.close()

#similar to what I want to do: 
#https://stackoverflow.com/questions/42982962/renaming-name-id-in-gffile


#end of script 
