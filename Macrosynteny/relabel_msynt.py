#! /usr/bin/env python3

import glob
import sys

try:
    _, labels, prefix_length = sys.argv
except ValueError:
    print("""usage: relabel_msynt.py mylabels.labels prefixlength
    globs all *.msynt files in the $CWD, "relab" infix added before msynt suffix
    prefixlength the length of the suffix in the chromosome name. should match with the one in msyntfile acessions""")
    sys.exit()

prefix_length = int(prefix_length)

label_dict = {}
with open(labels, 'r') as f:
    for line in f:
        newname, oldname, _ = line.rstrip().split()
        prefix = newname[0:prefix_length].upper()
        if label_dict.get(prefix) is None:
            label_dict[prefix] = {oldname: newname}
        else:
            label_dict[prefix][oldname] = newname

file_ls = glob.glob("*.msynt")
for msynt_file in file_ls:
    outfile = msynt_file.replace('msynt', 'relab.msynt')
    with open(msynt_file, 'r') as f, open(outfile, 'w') as g:
        for line in f:
            oldchr, acc, pos, ortho, para = line.rstrip().split()
            prefix = acc.split('_')[0]
            newchr = label_dict[prefix][oldchr]
            outstr = '\t'.join((newchr, acc, pos, ortho, para))+'\n'
            g.write(outstr)
