#!/usr/bin/env python3

import sys


if len(sys.argv) != 3:
    print("""Usage:
    cleanup_clus.py acc_list myclusfile.clus > cleaned_up.clus

    acc_list: list of accessions, newline separated
    myclusfile.clus: the clus file u want to look into

    The script deletes from all the provided OGs the accessions that are absent from the provided list
    In case of proteins that are in the list but not in the clus, a singleton OG is added
    Orthofinders pHOGs are splitted OGs, so post-pHOG singletons are in the output file (e.g. N0)
    However, pre-pHOG singletons are not included, because not informative orthology wise.
""")
    sys.exit()

_, acc_file, clus_file = sys.argv

acc_set = set()
with open(acc_file, 'r') as f:
    for acc in f:
        acc = acc.rstrip()
        acc_set.add(acc)

non_singletons = set()
with open(clus_file, 'r') as f:
    for line in f:
        og, _, *acc_ls = line.rstrip().split('\t')
        acc_ls_filt = [acc for acc in acc_ls if acc in acc_set]
        non_singletons.update(acc_ls_filt)
        acc_str = '\t'.join(acc_ls_filt)
        sys.stdout.write(f'{og}\t{len(acc_ls_filt)}\t{acc_str}\n')

singletons = acc_set - non_singletons

#Orthofinders pHOGS (even root ones) do not include singletons
#Could map to OG name (pre-pHOG), but since those are also singletons, it doesn't really matter
x = 0
for acc in singletons:
    x += 1
    sys.stdout.write(f'singleton{x:07}\t1\t{acc}\n')



