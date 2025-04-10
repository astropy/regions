# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pathlib import Path

import click
import pyregion
from astropy import log

from regions import read_ds9

log.setLevel('DEBUG')

TEST_FILE_DIR = Path('../regions/io/ds9/tests/data')


@click.group()
def cli():
    """
    Region parser debugging tool.
    """
    pass


@cli.command('list-files')
def list_files():
    print('Available files')
    for ffile in TEST_FILE_DIR.glob('*.reg'):
        print(ffile.parts[-1])


@cli.command('parse')
@click.option('--interactive', is_flag=True, default=False)
@click.option('--parser', default='regions')
@click.option('--errors', default='strict')
@click.argument('filename')
def parse(filename, interactive, parser, errors):
    readname = TEST_FILE_DIR / filename
    print(f'Reading {readname}')
    print(f'Using parser {parser}')
    if parser == 'regions':
        regions = read_ds9(str(readname), errors=errors)
    elif parser == 'pyregion':
        regions = pyregion.open(str(readname))
    print(regions)
    if interactive:
        import IPython
        IPython.embed()


if __name__ == '__main__':
    cli()
