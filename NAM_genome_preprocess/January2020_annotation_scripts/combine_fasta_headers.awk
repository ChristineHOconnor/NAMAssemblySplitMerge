#!/usr/bin/awk -f

/^>/    { header = $0 }
!/^>/   { sequence[header] = sequence[header] $0 }

END {
    for (head in sequence) {
        printf("%s\n%s\n", head, sequence[head])
    }
}

#script that will combine CDS sequences from the same transcript under 1 header
#Why is this necessary? I used bedtools getfasta to get a CDS fasta file. However, I had to use each CDS annotation (every exon etc) from the gff file so there were multiple entries for the same gene transcript.
#To run: awk -f ./script.awk file.fa > output_cds.fasta
#script copied from: https://unix.stackexchange.com/questions/339432/combine-multi-fasta-sequence
#Also: add sed command to remove strand information (-) or (+) from fasta header 
