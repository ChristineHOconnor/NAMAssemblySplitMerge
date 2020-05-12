#These are some commands I ran to test how we could possibly find all genes that map 1:1 across 26 different genomes 
#It will be a lot more complicated to expand this to 26 genomes, but I think join should work and, if you automatically generate the join command it will be more straight-forward than trying to write a python script 

join -1 1 -2 1 -a 1 -a 2 <(sort -k 1,1 test_one_to_one_B_M.txt) <(sort -k 1,1 test_one_to_one_B_P.txt) | tr " " "\t" | join -1 1 -2 1 -a 1 -a 2 <(sort -k 1,1 -) <(sort -k 1,1 test_one_to_one_B_W.txt) | tr " " "\t" > test_one_to_one_B_M_P_W.txt

#There are genes that link between B-M-P but don't show up in the B-P one to one comarison file. Most of these appear to be because the genes map as split/merge
