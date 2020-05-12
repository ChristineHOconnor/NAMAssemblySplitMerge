#! /usr/bin/perl -w

# pull_out_caononical_transcript.pl
# 28 April 2020
# Candy Hirsch

use strict;
use Getopt::Std;

my $usage = "\n$0 -i input_transcript_fasta_file -l input_list_of_canonical_transcripts -o output_longest_transcript_file\n";

our ($opt_i, $opt_l, $opt_o, $opt_h);
getopts("i:l:k:o:h") or die "$usage";

if ( (!(defined $opt_i)) || (!(defined $opt_l)) || (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
  exit;
}

open (my $FastaIn_fh, '<', $opt_i) || die "\nCannot open $opt_i\n\n";
open (my $list_fh, '<', $opt_l) || die "\nCannot open $opt_l\n\n";
open (my $FastaOut_fh, '>', $opt_o) || die "\nCannot open $opt_o\n\n";

my %canonical;
while (my $line = <$list_fh>) {
  chomp $line;
  $canonical{$line} = 0;
}

#Pull out information for first sequence
my $seq_name = <$FastaIn_fh>;
chomp $seq_name;

#Read through the fasta file
my $seq;
while (my $line = <$FastaIn_fh>) {
  chomp $line;
  if ($line !~ /^>/) {
    $seq .= "$line\n";
  }
  
  if ($line =~ /^>/) {


    if ($seq_name =~ /^>(\S+)/) {
      
      #Finish the previous sequence
      if (defined $canonical{$1}) {
	++$canonical{$1};
	print $FastaOut_fh "$seq_name\n$seq";
      }
    }

    undef $seq_name;
    undef $seq;
    $seq_name = $line;
  } 
}

#Finish the last sequence
if ($seq_name =~ /^>(\S+)/) {
  if (defined $canonical{$1}) {
    ++$canonical{$1};
    print $FastaOut_fh "$seq_name\n$seq";
  }
}

foreach my $keys (keys %canonical) {
  if ($canonical{$keys} != 1) {
    print "$keys\t$canonical{$keys}\n";
  }
}
