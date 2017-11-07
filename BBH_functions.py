#!/usr/bin/python
"""Functions for running comparing two proteomes and returning
Best Bilateral Hits using BLAST."""
import argparse
from os.path import isfile
from os.path import isdir
from os.path import abspath
from os import makedirs
import subprocess
import re
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

def getparms():
    """Retrieve arguments entered by the user."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-db1",
                        help="your input protein sequence in fasta format")
    parser.add_argument("-db2",
                        help="your input BLAST database")
    parser.add_argument("-out",
                        help="the name of the directory where you wish to store output")
    parser.add_argument("--evalue",
                        default=1E-30,
                        help="The minimum acceptable BLAST evalue")
    args = parser.parse_args()
    return args


def makeblastdb(db):
    """Simple function to make protein BLAST db from faa file.
       Input: path to a faa file.
       Output: no return value, just runs makeblastdb command."""
    #Does this file exist already?
    if isfile(db + ".psq"):
        pass
    #If it doesn't make us a protein blast DB.
    else:
        make2ndBLASTdbcmd = "makeblastdb -in " + db + " -input_type fasta -dbtype prot"
        subprocess.call(make2ndBLASTdbcmd, shell=True)

def getbbh(dict1, dict2):
    """Compare BLAST hits from db1->db2 vs. db2->db1;
       return Best Bilateral hits.
       Input:
       Dictionary containing hits from db1->db2 BLAST
       Dictionary containing hits from db2->db1 BLAST
       Returns:
       Dictionary with keys = db1 proteins,
       values = corresponding db2 proteins."""
    #pull out IDS that are in values of one and keys of another.
    #there may be a more elegant way to do this.
    BBH_dict = {}
    for k in dict1.keys():
        # v is the ID of the best hit of k from db1 into db2
        v = dict1[k]["hit ID"]
        # is k the ID of the best hit of v from db2 into db1???
        if  dict2[v]["hit ID"] == k:
            BBH_dict[k] = dict1[k]
    return BBH_dict

def printbbhs(iddict):
    """Simple function to print the results in tab deliminated form."""
    for k in iddict.keys():
        ID1, ID2 = k, iddict[k]["hit ID"]
        print(ID1 + "\t" + ID2)

def getmappings(mapfile):
    """I don't know what this is for...."""
    rhandle = open(mapfile, "r")
    mappings = {}
    for line in rhandle:
        ids = line.split("\t")
        mappings[ids[0]] = ids[1]
    return mappings

def blastindex(db1, db2, outputdir, evalue=1E-30):
    """Run blast search and parse the results.
       Input:
       db1 and db2 are BLAST databases.
       outputdir is a place to put the BLAST results.
       Returns:
       A dictionary where keys are the query proteins,
       values are the corresponding hits."""
    blastoutputdir = outputdir + "/BLAST_data" #Make dir
    if not isdir(blastoutputdir):
        makedirs(abspath(blastoutputdir))
    outputID = "/" + re.split('[\\\/.]+', db1)[-2] + "_vs_" + re.split('[\\\/.]+', db2)[-2] + ".XML"
    print(outputID)
    blastoutputfile1 = abspath(blastoutputdir + outputID)
    #Make blastdbs from fasta:
    makeblastdb(db1)
    makeblastdb(db2)
    # run the BLASTs; 1 into 2 then 2 into 1
    cline = NcbiblastpCommandline(query=db1,
                                  db=db2,
    							  evalue=evalue,
    							  outfmt=5,
    							  out=blastoutputfile1)
    cline()
    #Parse the BLAST
    result_handle = open(blastoutputfile1, "r")
    blastrecs = NCBIXML.parse(result_handle)
    iddict1 = {}
    for B in blastrecs:
        if B.alignments:
            ali = B.alignments[0]
            qID = B.query
            hID = ali.hit_def
            iddict1[qID] = {"hit ID" : hID, "alignmnent" : ali}
    return iddict1

def diamondindex(db1, db2, outputdir, evalue=1E-30):
    """Run DIAMOND search and parse the results.
       Input:
       db1 and db2 are BLAST databases.
       outputdir is a place to put the BLAST results.
       Returns:
       A dictionary where keys are the query proteins,
       values are the corresponding hits."""
    print("db1 is....  " + db1)
    print("db2 is....  " + db2)
    #Make dir to store results
    blastoutputdir = outputdir + "/BLAST_data"
    if not isdir(blastoutputdir):
        makedirs(abspath(blastoutputdir))
    #the weird re split stuff is so I can use the file names of DBs,
    #not paths to make the output file name.
    species1 = re.split('[\\\/.]+', db1)[-2]
    species2 = re.split('[\\\/.]+', db2)[-2]
    outputID = "/" + species1 + "_vs_" + species2 + ".XML"
    blastoutputfile = abspath(blastoutputdir + outputID)
    #Make blastdbs from fasta:
    makediamond = "diamond  makedb --in " + db2 + " --db " + db2
    print(makediamond)
    subprocess.call(makediamond, shell=True)
    rundiamond = ("diamond blastp --db " + db2 +
                               " --out " + blastoutputfile +
                               " --query " + db1 +
                               " --outfmt 5 -e 1e-3 --quiet")
    print(rundiamond)
    subprocess.call(rundiamond)
    #Parse the BLAST
    result_handle = open(blastoutputfile, "r")
    blastrecs = NCBIXML.parse(result_handle)
    iddict1 = {}
    #consider making the iddict1 an ordered dict to ensure correctness in synteny assessment.
    for B in blastrecs:
        if B.alignments:
            ali = B.alignments[0]
            qID = B.query
            hID = ali.hit_def
            #THis is only the score of the top HSP, not for the whole thing
            score = ali.hsps[0].score
            iddict1[qID] = {"hit ID " + species2     : hID,
                               "score "  + species2     : score}
    return iddict1
