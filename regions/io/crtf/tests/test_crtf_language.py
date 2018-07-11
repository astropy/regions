
import distutils.version as vers
import pytest


import astropy.version as astrov
from astropy.utils.data import get_pkg_data_filename


from ..read import CRTFParser, read_crtf
from ..write import crtf_objects_to_string
from ..core import CRTFRegionParserError

_ASTROPY_MINVERSION = vers.LooseVersion('1.1')
_ASTROPY_VERSION = vers.LooseVersion(astrov.version)

implemented_region_types = ('ellipse', 'circle', 'rectangle', 'poly', 'point', 'text', 'symbol')


def test_global_parser():
    """
    Checks that the global_parser does what's expected.
    """
    global_test_str = str("global coord=B1950_VLA, frame=BARY, corr=[I, Q], color=blue")
    global_parser = CRTFParser(global_test_str)
    assert dict(global_parser.global_meta) == {'coord': 'B1950_VLA', 'frame': 'BARY',
                                               'corr': ['I', 'Q'], 'color': 'blue'}


def test_valid_crtf_line():
    """
    Checks whether a the line is valid CRTF format.
    """
    line_str = 'coord=B1950_VLA, frame=BARY, corr=[I, Q], color=blue'

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(line_str)

    assert 'Not a valid CRTF line:' in str(err)


def test_valid_region_type():
    """
    Checks whether the region type is valid in CRTF format
    """
    reg_str = 'hyperbola[[18h12m24s, -23d11m00s], 2.3arcsec]'

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(reg_str)

    assert "Not a valid CRTF Region type: 'hyperbola'" in str(err)


def test_valid_global_meta_key():
    """
    Checks whether the global key is valid or not.
    """

    global_test_str = str("global label=B1950_VLA, frame=BARY, corr=[I, Q], color=blue")

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(global_test_str)

    assert "'label' is not a valid global meta key" in str(err)


def test_valid_meta_key():
    """
    Checks whether the key is valid or not.
    """

    meta_test_str = str("annulus[[17h51m03.2s, -45d17m50s], [0.10deg, 4.12deg]], hello='My label here'")

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(meta_test_str)

    assert "'hello' is not a valid meta key" in str(err)


def test_valid_region_syntax():
    """
    Checks whether the region has valid parameters
    """

    reg_str1 = "circle[[18h12m24s, -23d11m00s], [2.3arcsec,4.5arcsec]"

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(reg_str1)

    assert "Not in proper format: ('2.3arcsec', '4.5arcsec') should be a single length" in str(err)

    reg_str2 = "symbol[[32.1423deg, 12.1412deg], 12deg], linewidth=2, coord=J2000, symsize=2"

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(reg_str2)

    assert "Not in proper format: '12deg' should be a symbol" in str(err)

    reg_str3 = "circle[[18h12m24s, -23d11m00s]"

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(reg_str3)

    assert "Does not contain expected number of parameters for the region 'circle'" in str(err)

    reg_str4 = "poly[[1, 2], [3,4], [5, 6]]"

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(reg_str4)

    assert "polygon should have > 4 coordinates" in str(err)

    reg_str5 = "poly[[1, 2], [3,4], [5, 6],[1,6]]"

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(reg_str5)

    assert "In polygon, the last and first coordinates should be same" in str(err)

    reg_str6 = "rotbox[[12h01m34.1s, 12d23m33s], [3arcmin,], 12deg]"

    with pytest.raises(CRTFRegionParserError) as err:
        CRTFParser(reg_str6)

    assert "('3arcmin', '') should be a pair of length" in str(err)


@pytest.mark.parametrize('filename', ['data/CRTFgeneral.crtf', 'data/CRTFgeneraloutput.crtf'])
def test_file_crtf(filename):

    filename = get_pkg_data_filename(filename)
    regs = read_crtf(filename, 'warn')
    actual_output = crtf_objects_to_string(regs, 'fk4', '.3f').strip()

    with open(get_pkg_data_filename('data/CRTFgeneraloutput.crtf')) as f:
        ref_output = f.read().strip()

    # since metadata is not required to preserve order, we have to do a more
    # complex comparison
    desired_lines = [set(line.split(",")) for line in ref_output.split("\n")]
    actual_lines = [set(line.split(",")) for line in actual_output.split("\n")]
    for split_line in actual_lines:
        assert split_line in desired_lines

    for split_line in desired_lines:
        assert split_line in actual_lines
