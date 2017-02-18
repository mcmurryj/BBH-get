###BBH getter
#inputs: 2x genomes in fasta format
#outputs: a table of
#Get parameters
import argparse
import os
from Bio.Blast.Applications import NcbiblastpCommandline

def get_parms():
    parser = argparse.ArgumentParser()
    parser.add_argument("-db1",
                        help = "your input protein sequence in fasta format")
    parser.add_argument("-db2",
                        help = "your input BLAST database")
    parser.add_argument("-out",
                        help = "the name of the directory where you wish to store output")
    parser.add_argument("--evalue",
                        default = 1E-30,
                        help = "The minimum acceptable BLAST evalue")
    args = parser.parse_args()
    return(args)
    evalue     = args.evalue
    db1        = os.path.abspath(args.db1)
    db2        = os.path.abspath(args.db2)
    output_dir = os.path.abspath(args.out)

def get_BBH :
    #Make dir to store results
    blastoutputdir = output_dir + "/BLAST_data"
    os.makedirs(os.path.abspath(blastoutputdir))
    blastoutputfile1 = os.path.abspath(blastoutputdir + "/blast_output1.XML")
    blastoutputfile2 = os.path.abspath(blastoutputdir + "/blast_output2.XML")
    #Make blastdbs from fasta:
    if exists (args.db1 + "psq") and exists (args.db2 + "psq") :            #pseudocode must fix up like real
        pass
    else:
        make2ndBLASTdbcmd = "makeblastdb -in " + args.db1 + " -input_type fasta -dbtype prot"
        make2ndBLASTdbcmd = "makeblastdb -in " + args.db2 + " -input_type fasta -dbtype prot"
    # run the BLASTs; 1 into 2 then 2 into 1
    cline = NcbiblastpCommandline(query  = db1,
                                  db     = db2,
    							  evalue = evalue,
    							  outfmt = 5,
    							  out    = blastoutputfile1)
    cline()
    cline = NcbiblastpCommandline(query  = db2,
                                  db     = db1,
    							  evalue = evalue,
    							  outfmt = 5,
    							  out    = blastoutputfile2)
    cline()
    #Parse the BLAST
    result_handle        = open( blastoutputfile1, "r")
    BLASTrecs            = NCBIXML.parse(result_handle)
    ID_dict1 = {}
    for B in BLASTrecs :
        qID          = B.query_ID
        hID          = B.alignments[0].hit_ID
        ID_dict1[sorted([hID, qID])] = alignments[0].hsps[0].score

    ID_dict2 = {}
    for B in BLASTrecs :
        qID          = B.query_ID
        hID          = B.alignments[0].hit_ID
        ID_dict1[sorted([hID, qID])] = alignments[0]

    #pull out IDS that are in values of one and keys of another.
    oneTOtwo = set(ID_dict1.keys())
    twoTOone = set(ID_dict2.keys())
    BBHs     = list(oneTOtwo, twoTOone)
    return(ID_dict_1[BBHs])
