"""Compare DS9 parsing of the astropy regions package to pyregion

This scripts compares the DS9 parsing of the astropy regions package to
pyregion in two regards.

* speed : the time to parse a file given ``REPETITIONS`` repetitions

* completeness : the ratio of regions in a file (estimate by the number of '('
  in the file) that survives the parsing process

The test files are the ones used since PyAstro16, see
https://zenodo.org/record/5679://zenodo.org/record/56793 

"""
from pathlib import Path
from astropy.table import Table
from regions import read_ds9, write_ds9, DS9RegionParserError
import pyregion
import timeit
import numpy as np

TEST_FILE_DIR = Path('../regions/io/tests/data')
REPETITIONS = 1

results = list()

for filename in TEST_FILE_DIR.glob('*.reg'):
    print('\n\n{}'.format(filename))

    # estimate total number of regions
    n_regions = 0
    with open(str(filename)) as origin_file:
        for line in origin_file:
            if '(' in line:
                n_regions += 1

    # regions
    try:
        region_regions = read_ds9(str(filename))
    except DS9RegionParserError:
        time_regions = -1
        compl_regions = 0
    else:
        time_regions = timeit.timeit(
            'read_ds9(str(filename))',
            setup='from regions import read_ds9',
            globals=globals(),
            number=REPETITIONS) / REPETITIONS

        compl_regions = np.divide(len(region_regions), n_regions)

    # pyregion
    try:
        pyregion_regions = pyregion.open(str(filename))
    except ValueError:
        time_pyregion = -1
        compl_pyregion = 0
    else:
        time_pyregion = timeit.timeit(
            'pyregion.open(str(filename))',
            setup='import pyregion',
            globals=globals(),
            number=REPETITIONS) / REPETITIONS

        compl_pyregion = np.divide(len(pyregion_regions), n_regions)

    # collect results
    results.append(dict(
        filename=filename.parts[-1],
        time_regions=time_regions,
        time_pyregion=time_pyregion,
        compl_regions=compl_regions,
        compl_pyregion=compl_pyregion,
    ))


result_table = Table()
result_table['filename'] = [_['filename'] for _ in results]
result_table['time_regions'] = [_['time_regions'] for _ in results]
result_table['time_regions'].format = '.3f'
result_table['time_pyregion'] = [_['time_pyregion'] for _ in results]
result_table['time_pyregion'].format = '.3f'
result_table['compl_pyregion'] = [_['compl_pyregion'] for _ in results]
result_table['compl_pyregion'].format = '0.1%'
result_table['compl_regions'] = [_['compl_regions'] for _ in results]
result_table['compl_regions'].format = '0.1%'
result_table.meta['repetitions for timing'] = REPETITIONS
result_table.pprint()
