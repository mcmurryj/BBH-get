#!/usr/bin/python
from BBH_functions import get_parms
args                 = get_parms()
db1, db2, output_dir, evalue = args.db1, args.db2, args.out, args.evalue

from BBH_functions import BLAST_index
dict1 = BLAST_index(db1, db2, output_dir, evalue)
dict2 = BLAST_index(db2, db1, output_dir, evalue)

from BBH_functions import get_BBH
BBH_dict = get_BBH(dict1, dict2)

from BBH_functions import print_BBHs
print_BBHs(BBH_dict)
