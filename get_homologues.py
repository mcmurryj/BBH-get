#!/usr/bin/python
parser = argparse.ArgumentParser()
parser.add_argument("-hs",
                    help = "Subunits of human")
parser.add_argument("-sc",
                    help = "Subunits of yeast")
parser.add_argument("-db",
                    help = "your input BLAST database")
parser.add_argument("-out",
                    help = "the name of the directory where you wish to store output")
parser.add_argument("-map",
                    help = "TSV file mapping sc to hs subunits")
parser.add_argument("--evalue",
                    default = 1E-30,
                    help = "The minimum acceptable BLAST evalue")
args = parser.parse_args()

from BBH_functions import BLAST_index
#keys are hs subunits; values are info on hits in other genome.
dict_hs = BLAST_index(args.hs, args.db, args.output_dir, args.evalue)
# keys are sc subunits; values are info on hits in other genome.
dict_sc = BLAST_index(args.sc, args.db, args.output_dir, args.evalue)

#A dict with the
from BBH_functions import get_mappings
mapping_dict = get_mappings(args.map)

for k in dict_sc.keys():
    sc_hit_id = dict_sc[k]
    j         = mapping_dict[k]
    hs_hit_id = dict_hs[j]
    if sc_hit_id == j_hit_id :
        print "WOOO"

from BBH_functions import print_BBHs
print_BBHs(BBH_dict)
