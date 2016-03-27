from regions.io import read_ds9
import glob
import time

l = glob.glob('data/*.reg')

from regions.io.read_ds9 import ds9_parser, line_parser, region_list_to_objects, type_parser, coordinate, angular_length_quantity
# %lprun -f type_parser -f ds9_parser -f line_parser -f region_list_to_objects -f coordinate -f angular_length_quantity read_all_test_files()


def read_all_test_files():
    for f in l:
        print(f)
        if 'strip' in f or 'comment' in f:
            continue
           
        t0 = time.time()
        read_ds9(str(f))
        print(time.time()-t0)
