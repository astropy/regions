import pytest

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename

from ..core import CRTFRegionParserError
from ..read import CRTFParser, read_crtf
from ..write import crtf_objects_to_string
from ....shapes.circle import CircleSkyRegion
from ....shapes.ellipse import EllipseSkyRegion


implemented_region_types = ('ellipse', 'circle', 'rectangle', 'poly', 'point',
                            'text', 'symbol')


def test_global_parser():
    """
    Checks that the global_parser does what's expected.
    """
    global_test_str = ("global coord=B1950_VLA, frame=BARY, corr=[I, Q], "
                       "color=blue")
    global_parser = CRTFParser(global_test_str)
    expected = {'coord': 'B1950_VLA', 'frame': 'BARY',
                'corr': ['I', 'Q'], 'color': 'blue'}
    assert dict(global_parser.global_meta) == expected


def test_valid_crtf_line():
    """
    Checks whether a the line is valid CRTF format.
    """
    line_str = 'coord=B1950_VLA, frame=BARY, corr=[I, Q], color=blue'

    with pytest.raises(CRTFRegionParserError) as excinfo:
        CRTFParser(line_str)

    assert 'Not a valid CRTF line:' in str(excinfo.value)


def test_valid_region_type():
    """
    Checks whether the region type is valid in CRTF format
    """
    reg_str = 'hyperbola[[18h12m24s, -23d11m00s], 2.3arcsec]'

    with pytest.raises(CRTFRegionParserError) as excinfo:
        CRTFParser(reg_str)

    assert "Not a valid CRTF Region type: 'hyperbola'" in str(excinfo.value)


def test_valid_global_meta_key():
    """
    Checks whether the global key is valid or not.
    """
    global_test_str = ("global label=B1950_VLA, frame=BARY, corr=[I, Q], "
                       "color=blue")
    with pytest.raises(CRTFRegionParserError) as excinfo:
        CRTFParser(global_test_str)
    assert "'label' is not a valid global meta key" in str(excinfo.value)


def test_valid_meta_key():
    """
    Checks whether the key is valid or not.
    """
    meta_test_str = ("annulus[[17h51m03.2s, -45d17m50s], [0.10deg, 4.12deg]], "
                     "hello='My label here'")
    with pytest.raises(CRTFRegionParserError) as excinfo:
        CRTFParser(meta_test_str)
    assert "'hello' is not a valid meta key" in str(excinfo.value)


def test_valid_region_syntax():
    """
    Checks whether the region has valid parameters
    """
    reg_str1 = "circle[[18h12m24s, -23d11m00s], [2.3arcsec,4.5arcsec]"
    with pytest.raises(CRTFRegionParserError) as excinfo:
        CRTFParser(reg_str1)
    assert ("Not in proper format: ('2.3arcsec', '4.5arcsec') should be "
            "a single length" in str(excinfo.value))

    reg_str2 = ("symbol[[32.1423deg, 12.1412deg], 12deg], linewidth=2, "
                "coord=J2000, symsize=2")
    with pytest.raises(CRTFRegionParserError) as excinfo:
        CRTFParser(reg_str2)
    assert ("Not in proper format: '12deg' should be a symbol" in
            str(excinfo.value))

    reg_str3 = "circle[[18h12m24s, -23d11m00s]"
    with pytest.raises(CRTFRegionParserError) as excinfo:
        CRTFParser(reg_str3)
    assert ("Does not contain expected number of parameters for the region "
            "'circle'" in str(excinfo.value))

    reg_str4 = "poly[[1, 2], [4, 5]]"
    with pytest.raises(CRTFRegionParserError) as excinfo:
        CRTFParser(reg_str4)
    assert "polygon should have >= 3 coordinates" in str(excinfo.value)

    reg_str6 = "rotbox[[12h01m34.1s, 12d23m33s], [3arcmin,], 12deg]"
    with pytest.raises(CRTFRegionParserError) as excinfo:
        CRTFParser(reg_str6)
    assert "('3arcmin', '') should be a pair of length" in str(excinfo.value)


def test_issue_312_regression():
    """
    Make sure there is no trailing comma when writing a CRTF string with
    no metadata.
    """
    reg = EllipseSkyRegion(center=SkyCoord(279.174990 * u.deg,
                                           -21.257123 * u.deg, frame='fk5'),
                           width=0.001571 * u.deg, height=0.001973 * u.deg,
                           angle=111.273322 * u.deg)
    crtfstr = crtf_objects_to_string([reg], 'fk5', '.6f', 'deg')
    assert crtfstr.strip()[-1] != ','


@pytest.mark.parametrize('filename', ['data/CRTFgeneral.crtf',
                                      'data/CRTFgeneraloutput.crtf'])
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


def test_crtf_header():
    """
    Regression test for issue #239.
    """
    crtf_str = ('#CRTFv0 CASA Region Text Format version 0\n'
                'circle[[42deg, 43deg], 3deg], coord=J2000, color=green')
    parser = CRTFParser(crtf_str)
    reg = parser.shapes.to_regions()[0]
    assert isinstance(reg, CircleSkyRegion)
    assert reg.center.ra.value == 42.0
    assert reg.center.ra.unit == 'deg'
    assert reg.center.dec.value == 43.0
    assert reg.center.dec.unit == 'deg'
    assert reg.radius.value == 3.0
    assert reg.radius.unit == 'deg'


def test_space_after_regname():
    """
    Regression test for #271: space is allowed
    """
    reg_str = 'circle [[42deg, 43deg], 3deg], coord=J2000, color=green'
    parser = CRTFParser(reg_str)
    reg = parser.shapes.to_regions()[0]
    assert isinstance(reg, CircleSkyRegion)


def test_no_comma_after_region():
    reg_str = 'circle [[42deg, 43deg], 3deg] coord=J2000, color=green'
    parser = CRTFParser(reg_str)
    reg = parser.shapes.to_regions()[0]
    assert isinstance(reg, CircleSkyRegion)
    assert reg.center.ra.value == 42.0
    assert reg.center.ra.unit == 'deg'
    assert reg.center.dec.value == 43.0
    assert reg.center.dec.unit == 'deg'
    assert reg.radius.value == 3.0
    assert reg.radius.unit == 'deg'
