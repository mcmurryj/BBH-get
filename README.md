get_BBH is a tool for retrieving best bilateral blast hits between two fasta databases.

USAGE:
get_BBH -db1 your_first_faa_file.faa -db2 your_2nd_FAA_file.faa -out /some/path/youroutputdir --evalue 1E-50

Output:
right now, prints the IDs of best bilateral BLAST hit pairs.  Pretty messy.  

TODO:
Finish the get_homologues script, which will retrieve sequences with the same best hit from 2x reference sets.  For assembling MSA of large complex type stuff.
