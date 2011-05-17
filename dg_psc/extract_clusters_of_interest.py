#!/usr/bin/env python
"""Print info on and rearrange clusters around proteins of interest.

This automates boring by-hand examination of output clustered HTML files to help
facilitate rapid cycles of examining different clustering parameters.
"""
from __future__ import with_statement
import sys
import os
import glob
import shutil

interest_items = dict(
  IPR000953 = [('Cbx2', 'P30658', True),
               ('Cbx8', 'Q9QXV1', False),
               ('Cbx4', 'O55187', False),
               ('Cbx6', 'Q9DBY5', False),
               ('Cbx7', 'Q8VDS3', False),
               ('Pc', 'P26017', True)],
  IPR001841 = [('Psc', 'P35820', True),
               ('Suz2', 'P25172', False),
               ('Bmi1', 'P25916', True),
               ('Mel18', 'P23798', False),
               ('l_3_73Ah', 'Q9VV77', False)]
)

def main(ipr_number, work_dir):
    #out_dir = os.path.join(work_dir, "no_interactions_cluster")
    out_dir = work_dir
    new_files = glob.glob(os.path.join(work_dir, "%s-cluster*.html" %
        ipr_number))
    if len(new_files) > 0:
        #remove_old_files(out_dir, ipr_number)
        #for new_file in new_files:
        #    shutil.move(new_file, out_dir)
        examine_files(out_dir, interest_items[ipr_number], ipr_number)

def examine_files(base_dir, interest_items, ipr_number):
    all_files = glob.glob(os.path.join(base_dir, "%s-cluster*.html" %
        ipr_number))
    rename_list = []
    for cur_file in all_files:
        for short_name, uniprot_id, do_move in interest_items:
            with open(cur_file) as cur_handle:
                for line in cur_handle:
                    if line.find(uniprot_id) >= 0:
                        print short_name, os.path.split(cur_file)[-1]
                        if do_move:
                            new_name = os.path.join(base_dir,
                              "%s-cluster-%s.html" % (ipr_number, short_name))
                            rename_list.append((cur_file, new_name))
                        break
    for old_file, new_file in rename_list:
        shutil.move(old_file, new_file)

def remove_old_files(base_dir, ipr_number):
    for to_remove in glob.glob(os.path.join(base_dir, "%s-cluster*.html" %
        ipr_number)):
        os.remove(to_remove)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
