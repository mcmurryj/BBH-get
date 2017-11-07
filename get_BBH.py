#!/usr/bin/python
from BBH_functions import getparms
args                         = getparms()
db1, db2, output_dir, evalue = args.db1, args.db2, args.out, args.evalue

from BBH_functions import blastindex
dict1 = blastindex(db1, db2, output_dir, evalue)
dict2 = blastindex(db2, db1, output_dir, evalue)

from BBH_functions import getbbh
BBH_dict = getbbh(dict1, dict2)

from BBH_functions import printbbhs
printbbhs(BBH_dict)
