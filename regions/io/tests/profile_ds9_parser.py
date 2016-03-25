from regions.io import read_ds9
import glob

l = glob.glob('data/*.reg')

for f in l:
    print(f)
    if 'strip' in f or 'comment' in f:
        continue
       
    read_ds9(str(f))




