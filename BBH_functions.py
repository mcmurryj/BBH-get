#!/usr/bin/python
import argparse
import os
from Bio.Blast.Applications import NcbiblastpCommandline
import subprocess
from Bio.Blast import NCBIXML
import re

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
    if not os.path.isdir(blastoutputdir) :
        os.makedirs(os.path.abspath(blastoutputdir))
    outputID = "/" + re.split('[\\\/.]+', db1)[-2] + "_vs_" + re.split('[\\\/.]+', db2)[-2] + ".XML"
    print(outputID)
    blastoutputfile1 = os.path.abspath(blastoutputdir + outputID)
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

def get_BBH(dict1, dict2) :
    #pull out IDS that are in values of one and keys of another.
    #there may be a more elegant way to do this.
    BBH_dict = {}
    for k in dict1.keys() :
        # v is the ID of the best hit of k from db1 into db2
        v = dict1[k]["hit ID"]
        # is k the ID of the best hit of v from db2 into db1???
        if  dict2[v]["hit ID"] == k :
            BBH_dict[k] = dict1[k]
    return(BBH_dict)

def print_BBHs(ID_dict) :
    for k in ID_dict.keys() :
        ID1, ID2 = k, ID_dict[k]["hit ID"]
        print(ID1 + "\t" + ID2)

def get_mappings(mapfile) :
    rhandle  =  open(mapfile, "r")
    mappings =  {}
    for line in rhandle :
        ids  = line.split("\t")
        mappings[ids[0]] = ids[1]

    return(mappings)
