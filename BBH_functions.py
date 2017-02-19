#!/usr/bin/python
import argparse
import os
from Bio.Blast.Applications import NcbiblastpCommandline
import subprocess
from Bio.Blast import NCBIXML

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

def BLAST_index (db1, db2, output_dir, evalue = 1E-30) :
    #Make dir to store results
    blastoutputdir = output_dir + "/BLAST_data"
    os.makedirs(os.path.abspath(blastoutputdir))
    blastoutputfile1 = os.path.abspath(blastoutputdir + "/" +
                                       db1.split('.')[0] + "vs" +
                                       db2.split('.')[0] + ".XML")
    #Make blastdbs from fasta:
    if os.path.isfile(db2 + ".psq") :            #pseudocode must fix up like real
        pass
    else:
        make2ndBLASTdbcmd = "makeblastdb -in " + db2 + " -input_type fasta -dbtype prot"
        subprocess.call(make2ndBLASTdbcmd, shell = True)

    # run the BLASTs; 1 into 2 then 2 into 1
    cline = NcbiblastpCommandline(query  = db1,
                                  db     = db2,
    							  evalue = evalue,
    							  outfmt = 5,
    							  out    = blastoutputfile1)
    cline()
    #Parse the BLAST
    result_handle        = open( blastoutputfile1, "r")
    BLASTrecs            = NCBIXML.parse(result_handle)
    ID_dict1 = {}
    for B in BLASTrecs :
        if B.alignments :
            ali             = B.alignments[0]
            qID             = B.query
            hID             = ali.hit_def
            ID_dict1[qID]   = {"hit ID" : hID, "alignmnent" : ali}

    return(ID_dict1)

def get_BBH(dict1, dict2)
    #pull out IDS that are in values of one and keys of another.
    #there may be a more elegant way to do this.
    BBH_dict = {}
    for k in ID_dict1.keys() :
        # v is the ID of the best hit of k from db1 into db2
        v = ID_dict1[k]["hit ID"]
        # is k the ID of the best hit of v from db2 into db1???
        if  ID_dict2[v]["hit ID"] == k :
            BBH_dict[k] = ID_dict1[k]
    return(BBH_dict)

def print_BBHs(ID_dict) :
    for k in ID_dict.keys() :
        ID1, ID2 = k, ID_dict[k]["hit ID"]
        print(ID1 + "\t" + ID2)
