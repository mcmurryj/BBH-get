## Get_BBH

#### Description
Get_BBH is a tool for retrieving Best Bilateral BLAST hits (BBHs) between two fasta databases.

#### Usage
*get_BBH -db1 your_first_faa_file.faa -db2 your_2nd_FAA_file.faa -out /some/path/youroutputdir --evalue 1E-50*

#### Arguments
*db1* A proteome in multiple sequence fasta format.
*db2* A proteome in multiple sequence fasta format.
*out* A directory in which to store the BLAST output.
*evalue* Optional E value cutoff for BLAST search, default is 1E-30 which is stringent.

#### Output

Prints the IDs of best bilateral BLAST hit pairs to STDOUT.
