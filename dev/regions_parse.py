import click
import pyregion
from regions import read_ds9
from pathlib import Path

TEST_FILE_DIR = Path('../regions/io/tests/data')

@click.group()
def cli():
    """astropy.regions parser debugging tool."""
    pass

@cli.command('list-files')
def list_files():
    print("Available files")
    for ffile in TEST_FILE_DIR.glob('*.reg'):
        print(ffile.parts[-1])

@cli.command('parse')
@click.option('--interactive', is_flag=True, default=False)
@click.option('--parser', default='regions')
@click.argument('filename')
def parse(filename, interactive, parser):
    readname = TEST_FILE_DIR / filename
    print('Reading {}'.format(readname))
    print('Using parser {}'.format(parser))
    if parser == 'regions':
        regions = read_ds9(str(readname), errors='warn')
    elif parser == 'pyregion':
        regions = pyregion.open(str(readname))
    print(regions)
    if interactive:
        import IPython
        IPython.embed()

if __name__ == '__main__':
    cli()
