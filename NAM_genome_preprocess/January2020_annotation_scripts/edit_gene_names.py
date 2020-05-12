"This script will take a gff file and a file with current and new gene names and replace the old gene names in the gff file. A new gff file will be created and saved. Regular expressions are used to find and replace gene names. Originally I used python's replace function but that caused isssues becase searching for g7 would lead to replacements for g7 and also g70 etc." 

import sys
import getopt
import operator
import re

#usage etc setup taken from Meesh's script
def usage():
    print("""\n
        python3 edit_gene_name.py
            -g or gff_file gff file with original version of names
            -o or output_file  Name of new gff file that will be created 
            -k or name_key File with current and new gene names
        \n""")
try:
    opts, args = getopt.getopt(sys.argv[1:], "g:o:k:h", ["gff_file=",
                                                         "output_file=",
                                                         "name_key="
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
pattern_Zm = re.compile(r'(?<=gene:)(.*?)(?=;)|(?<=Name=)(.*?)(?=_)|(?<=transcript:)(.*?)(?=_)|(?<=CDS:)(.*?)(?=_)|(?<=transcript_id=)(.*?)(?=_)|(?<=exon_id=)(.*?)(?=_)|(?<=protein_id=)(.*?)(?=_)|(?<=transcript:)(.*?)(?=$)')
pattern_abi = re.compile(r'(?<=ID=)(.*?)(?=\.)|(?<=ID=)(.*?)(?=$)|(?<=Parent=)(.*?)(?=\.)|(?<=Parent=)(.*?)(?=$)')
#there's a problem with my pattern searches, so I'll add 1 at a time until I get an issue
#problem: I didn't escape '.' 
#new problem: Zm pattern's aren't being found 
#01/21/2020: I decided to create 2 pattern search's, 1 for Zm gene names and 1 for ab initio g* gene names
#this works on test data and on the full gff file
#I'm not sure why, but when I have all the regular expressions in 1 re.compile statement not everything is found 

with open(name_file) as gene_names,\
     open(gff_file) as f,\
     open(output_file, 'w') as out_f:
    # A dictionary mapping gene names to what they will be changed to
    rename_dict = dict(line.split() for line in gene_names)
    #print(rename_dict) #this is fine
    for line in f:
        if line.startswith("#") or 'assembly' in line:
            out_f.write(line)
        else: 
            # Search for the gene name
            #replace all instances of 'gene:' with nothing
            #new_line = line.replace("gene:","")
            #new_line2 = new_line.replace("transcript:", "")
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
#close the output file
out_f.close()

#try using sed to change all instances of 'gene:' and 'transcript:' to nothing
#subprocess.call()
#now replace all instances of 'gene:' and 'transcript:' with nothing
#Why? Because extra characters between ID= or Parent= and start of gene name will interfere with other steps in the pipeline
with open(output_file, "r") as sources: 
    lines = sources.readlines()
with open(output_file, "w") as sources:
    for line in lines: 
        sources.write(re.sub(r'=gene:|=transcript:', '=', line))
        #sources.write(re.sub(r'transcript:', '', line))

sources.close()

#similar to what I want to do: 
#https://stackoverflow.com/questions/42982962/renaming-name-id-in-gffile


#end of script 
