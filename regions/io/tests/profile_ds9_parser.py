from regions.io import read_ds9
import glob
import time

l = glob.glob('data/*.reg')


def read_all_test_files():
    for f in l:
        print(f)
        if 'strip' in f or 'comment' in f:
            continue

        t0 = time.time()
        read_ds9(str(f))
        print(time.time()-t0)
