#!/usr/bin/python
from BBH_getter import get_parms
args                 = get_parms()
db1, db2, output_dir, evalue = args.db1, args.db2, args.out, args.evalue

from BBH_getter import get_BBH
BBH_dict = get_BBH(db1, db2, output_dir, evalue)

from BBH_getter import print_BBHs
print_BBHs(BBH_dict)
