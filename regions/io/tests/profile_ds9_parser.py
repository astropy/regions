from regions.io import read_ds9
import glob

l = glob.glob('data/*.reg')

from regions.io.read_ds9 import ds9_parser, line_parser, region_list_to_objects
# %lprun -f ds9_parser -f line_parser -f region_list_to_objects read_all_test_files()


def read_all_test_files():
    for f in l:
        print(f)
        if 'strip' in f or 'comment' in f:
            continue
           
        read_ds9(str(f))
