# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the ds9 subpackage.
"""

import os

from astropy.coordinates import Angle, SkyCoord
from astropy.tests.helper import assert_quantity_allclose
import astropy.units as u
from astropy.utils.data import get_pkg_data_filenames
from astropy.utils.exceptions import AstropyUserWarning
from numpy.testing import assert_allclose
import pytest

from ....core import Regions
from ....shapes.circle import CircleSkyRegion
from ...._utils.optional_deps import HAS_MATPLOTLIB  # noqa


def test_read():
    # Check that all test files are readable
    filenames = get_pkg_data_filenames('data')
    with pytest.warns(AstropyUserWarning):
        for filename in filenames:
            Regions.read(filename, format='ds9')


def test_serialize():
    """
    Simple test for serialize.
    """
    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)
    expected = ('# Region file format: DS9 astropy/regions\nicrs\n'
                'circle(42.0000,43.0000,3.0000)\n')
    actual = region.serialize(format='ds9', precision=4)
    assert actual == expected


def test_parse():
    """
    Simple test for parse.
    """
    region_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
                  'circle(42.0000,43.0000,3.0000)\n')
    region = Regions.parse(region_str, format='ds9')[0]

    assert_allclose(region.center.ra.deg, 42)
    assert_allclose(region.center.dec.deg, 43)
    assert_quantity_allclose(region.radius, 3 * u.deg)


def test_exclude():
    """
    Simple parse test for an excluded region.
    """
    region_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
                  ' -circle(42.0000,43.0000,3.0000)\n')
    region = Regions.parse(region_str, format='ds9')[0]
    assert not region.meta['include']


def test_read_write(tmpdir):
    """
    Simple test for write and read.
    """
    center = SkyCoord(42, 43, unit='deg', frame='fk5')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)
    region.meta['text'] = 'ExampleText'

    filename = os.path.join(str(tmpdir), 'ds9.reg')
    region.write(filename, format='ds9')
    region2 = Regions.read(filename, format='ds9')[0]

    assert_quantity_allclose(region.center.ra, 42 * u.deg)
    assert_quantity_allclose(region.center.dec, 43 * u.deg)
    assert_quantity_allclose(region.radius, 3 * u.deg)
    assert 'text' in region2.meta
    assert region2.meta['text'] == 'ExampleText'


def test_invalid_region_warns():
    ds9_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
               'circle(42.0000,43.0000,3.0000)\ninvalidregion(blah)')

    with pytest.warns(AstropyUserWarning) as warn_results:
        regions = Regions.parse(ds9_str, format='ds9')
    assert len(regions) == 1
    assert len(warn_results) == 1


def test_global_parser():
    """
    Test parsing global metadata.
    """
    region_str = ('global color=green dashlist=8 3 width=1 '
                  'font="helvetica 10 normal roman" select=1 highlite=1 '
                  'dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
                  'icrs;circle(42.0000,43.0000,3.0000)')
    region = Regions.parse(region_str, format='ds9')[0]

    meta_expected = {'select': 1, 'highlite': 1, 'fixed': 0, 'edit': 1,
                     'move': 1, 'delete': 1, 'include': 1, 'source': 1}
    visual_expected = {'linewidth': 1, 'fontname': 'helvetica',
                       'fontsize': 10, 'fontweight': 'normal',
                       'fontstyle': 'normal', 'facecolor': 'green',
                       'edgecolor': 'green', 'default_style': 'ds9'}

    assert region.meta == meta_expected
    assert region.visual == visual_expected


def test_meta_color():
    region_str = ('# Region file format: DS9 astropy/regions\nicrs\n'
                  'circle(42.0000,43.0000,3.0000) # color=green\n'
                  'circle(43.0000,43.0000,3.0000) # color=orange\n')
    regions = Regions.parse(region_str, format='ds9')

    assert regions[0].visual['facecolor'] == 'green'
    assert regions[0].visual['edgecolor'] == 'green'
    assert regions[1].visual['facecolor'] == 'orange'
    assert regions[1].visual['edgecolor'] == 'orange'


def test_meta_color_override_global():
    """
    Color parsing test in the presence of a global color.
    """
    region_str = ('# Region file format: DS9 astropy/regions\n'
                  'global color=blue dashlist=8 3 width=1 '
                  'font="helvetica 10 normal roman" select=1 highlite=1 '
                  'dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
                  'icrs\ncircle(42.0000,43.0000,3.0000) # color=green\n'
                  'circle(42.0000,43.0000,5.0000) # color=orange\n'
                  'circle(42.0000,43.0000,7.0000)')
    regions = Regions.parse(region_str, format='ds9')
    assert regions[0].visual['facecolor'] == 'green'
    assert regions[1].visual['facecolor'] == 'orange'
    assert regions[2].visual['facecolor'] == 'blue'


def test_issue134_regression():
    region_str = 'galactic; circle(+0:14:26.064,+0:00:45.206,30.400")'
    region = Regions.parse(region_str, format='ds9')[0]
    assert_quantity_allclose(region.radius, 30.4 * u.arcsec)


def test_issue65_regression():
    region_str = 'J2000; circle 188.5557102 12.0314056 1" # color=red'
    region = Regions.parse(region_str, format='ds9')[0]
    assert_quantity_allclose(region.center.ra, 188.5557102 * u.deg)
    assert_quantity_allclose(region.center.dec, 12.0314056 * u.deg)
    assert_quantity_allclose(region.radius, 1 * u.arcsec)


def test_frame_uppercase():
    """
    Regression test for issue #236 (PR #237).
    """
    region_str = ('GALACTIC\ncircle(188.5557102,12.0314056,7.1245) # '
                  'text={This message has both a " and : in it} '
                  'textangle=30')
    region = Regions.parse(region_str, format='ds9')[0]
    assert_quantity_allclose(region.center.l, 188.5557102 * u.deg)
    assert_quantity_allclose(region.center.b, 12.0314056 * u.deg)
    assert_quantity_allclose(region.radius, 7.1245 * u.deg)


def test_pixel_angle():
    """
    Check whether angle in PixelRegions is a u.Quantity object.
    """
    region_str = 'image\nbox(1.5,2,2,1,0)'
    region = Regions.parse(region_str, format='ds9')[0]
    assert isinstance(region.angle, u.Quantity)
    assert_quantity_allclose(region.angle, 0. * u.deg)


def test_parse_formats():
    region1_str = 'image\ncircle(1.5, 2, 2)'
    region1 = Regions.parse(region1_str, format='ds9')[0]
    region2_str = 'image\ncircle(1.5i, 2i, 2i)'
    region2 = Regions.parse(region2_str, format='ds9')[0]
    assert region1 == region2

    region3_str = 'icrs\ncircle(1.5, 2, 2)'
    region3 = Regions.parse(region3_str, format='ds9')[0]
    region4_str = 'icrs\ncircle(1.5d, 2d, 2d)'
    region4 = Regions.parse(region4_str, format='ds9')[0]
    assert region3 == region4

    region5_str = 'icrs\ncircle(1.5r, 1r, 2r)'
    region5 = Regions.parse(region5_str, format='ds9')[0]
    assert_quantity_allclose(region5.center.ra, 1.5 * u.radian)
    assert_quantity_allclose(region5.center.dec, 1.0 * u.radian)
    assert_quantity_allclose(region5.radius, 2.0 * u.radian)

    region6_str = 'icrs\ncircle(1:20:30, 2:3:7, 2)'
    region6 = Regions.parse(region6_str, format='ds9')[0]
    assert_quantity_allclose(region6.center.ra, 20.124999999999996 * u.deg)
    assert_quantity_allclose(region6.center.dec, 2.051944444444444 * u.deg)
    assert_quantity_allclose(region6.radius, 2. * u.deg)

    region7_str = 'icrs\ncircle(1h20m30s, 2d3m7s, 2d)'
    region7 = Regions.parse(region7_str, format='ds9')[0]
    assert_quantity_allclose(region7.center.ra, 20.124999999999996 * u.deg)
    assert_quantity_allclose(region7.center.dec, 2.051944444444444 * u.deg)
    assert_quantity_allclose(region7.radius, 2. * u.deg)


def test_text_metadata():
    """
    Regression test for issue #233: make sure that text metadata is
    parsed and stored appropriately.
    """
    region_str = 'image\ncircle(1.5, 2, 2) # text={this_is_text}'
    regions = Regions.parse(region_str, format='ds9')
    assert len(regions) == 1
    assert regions[0].meta['text'] == 'this_is_text'


def test_angle_serialization():
    """
    Regression test for issue #223 to ensure Angle arcsec inputs are
    correctly converted to degrees.
    """
    region = CircleSkyRegion(SkyCoord(10, 20, unit='deg'), Angle(1, 'arcsec'))
    region_str = region.serialize(format='ds9', precision=6)

    expected = ('# Region file format: DS9 astropy/regions\nicrs\n'
                'circle(10.000000,20.000000,0.000278)\n')
    assert region_str == expected


def test_semicolon():
    """
    Regression test for issue #238 to allow semicolons in the text
    field.
    """
    region_str = ('galactic\n'
                  'circle(0.003, 0.1, 206.696") # text={S17; test} color=red\n'
                  'circle(0.1, -0.5, 360.148") # text={S19} color=red')

    regions = Regions.parse(region_str, format='ds9')
    assert len(regions) == 2

    regstr = ('galactic;'
              'circle(202.4,47.2,10.9) # text={test; test; A1} color=red;'
              'circle(103.8,40.6,3.6) # text="A, \'B\', ; {C}, D"'
              ' color=orange;'
              'circle(103.8,30.6,3.6) # text="A, \'B, ; C}, D" color=orange;'
              'circle(303.8,40.6,3.6) # text=\'A, "B", ; {C}, D\' color=cyan;'
              'circle(303.8,40.6,3.11);'
              'circle(303.8,30.6,3.6) # text=\'A, "B, ; {C, D\' color=cyan;'
              'circle(503.8,40.6,3.6) # text={A, "B", ; \'C\', D}'
              ' color=yellow\n'
              'circle(503.8,30.6,3.6) # text={A, "B, ; \'C, D} color=yellow;'
              '-ellipse(104.84,68.55,7.93",3.96",158.39) # dashlist=8 3'
              ' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1'
              ' source=1 text={Ellipse;Text} background color=#0ff width=1'
              ' font="helvetica 10 normal italic";'
              "ellipse(104.84,68.55,7.93',3.96',158.39) # dashlist=8 3"
              ' text="hello; world;";'
              'circle(1.5, 2, 2) # text={this_is_text};'
              'text(151.1,51.2) # text={text ; hi ; there;} textangle=45'
              ' color=orange font="helvetica 12 normal roman";'
              'text(151.1,51.2) # text  =  {text ; hi ; there;} textangle=45'
              ' color=orange font="helvetica 12 normal roman";'
              'circle(1.5, 2, 2) # text={text="hi;world";text;"tes\'t};'
              "circle(1.5, 2, 2) # text={text='hi;world';text;\"tes't};")

    regions = Regions.parse(regstr, format='ds9')
    assert len(regions) == 15


def test_parser_no_metadata():
    """
    Regression test for issue #259 to ensure regions without metadata
    are parsed correctly.
    """
    region1_str = 'galactic;circle(42,43,3)'
    region2_str = 'galactic;circle 42 43 3'
    region1 = Regions.parse(region1_str, format='ds9')[0]
    region2 = Regions.parse(region2_str, format='ds9')[0]

    assert isinstance(region1, CircleSkyRegion)
    assert isinstance(region2, CircleSkyRegion)
    assert_quantity_allclose(region1.radius, 3.0 * u.deg)
    assert_quantity_allclose(region2.radius, 3.0 * u.deg)


def test_spaces_metadata():
    """
    Regression test for spaces before or after "=" in region metadata.
    """
    regstr = ('galactic;'
              'circle(202.4,47.2,10.9) # text={a; test; test; A1} color=red;'
              'circle(202.4,47.2,10.9) # text ={a; test; test; A1} color= red;'
              'circle(202.4,47.2,10.9) # text= {a; test; test; A1} color =red;'
              'circle(202.4,47.2,10.9) # text = {a; test; test; A1}'
              ' color = red;')
    regions = Regions.parse(regstr, format='ds9')
    assert len(regions) == 4
    for i in range(1, len(regions)):
        assert regions[0].visual == regions[i].visual
        assert regions[0].meta == regions[i].meta

    regstr2 = ('galactic;'
               'text(151.1,51.2) # text={text ; hi ; there;} textangle= 45'
               ' color=orange font="helvetica 12 normal roman";'
               'text(151.1,51.2) # text  ={text ; hi ; there;} textangle =45'
               ' color=orange font="helvetica 12 normal roman";'
               'text(151.1,51.2) # text=  {text ; hi ; there;} textangle  =45'
               ' color   =   orange font="helvetica 12 normal roman";'
               'text(151.1,51.2) # text   =  {text ; hi ; there;} textangle=45'
               ' color   =orange font="helvetica 12 normal roman";')

    regions = Regions.parse(regstr2, format='ds9')
    assert len(regions) == 4
    for i in range(1, len(regions)):
        assert regions[0].visual == regions[i].visual
        assert regions[0].meta == regions[i].meta


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_point_boxcircle():
    import matplotlib.path as mpath

    region_str = ('# Region file format: DS9 astropy/regions\n'
                  'image\n'
                  'point(101.000000,101.000000) # color=red\n'
                  'point(101.000000,201.000000) # point=boxcircle color=blue\n'
                  'point(101.000000,301.000000)')
    regions = Regions.parse(region_str, format='ds9')

    assert isinstance(regions[0].as_artist().get_marker(), mpath.Path)
    assert regions[0].as_artist().get_markeredgecolor() == 'red'
    assert isinstance(regions[1].visual['marker'], mpath.Path)
    assert regions[1].as_artist().get_markeredgecolor() == 'blue'
    assert isinstance(regions[2].as_artist().get_marker(), mpath.Path)
    assert regions[2].as_artist().get_markeredgecolor() == '#00ff00'


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_compound_color():
    region_str = ('# Region file format: DS9 astropy/regions\n'
                  'image\n'
                  'annulus(651.0,301.0,60.0,90.0) # color=red')
    regions = Regions.parse(region_str, format='ds9')
    assert regions[0].as_artist().get_edgecolor() == (1., 0., 0., 1.)


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_default_mpl_kwargs():
    region_str = ('# Region file format: DS9 astropy/regions\n'
                  'image\n'
                  'circle(101.0,101.0,3.0) # color=red\n'
                  'circle(101.0,301.0,3.0)\n'
                  'point(101.0,101.0) # color=red\n'
                  'point(101.0,301.0)\n'
                  'text(101.0,101.0) # text={Text} color=red\n'
                  'text(101.0,101.0) # text={Text}')
    regions = Regions.parse(region_str, format='ds9')

    # Patch
    assert regions[0].visual['default_style'] == 'ds9'
    assert regions[0].as_artist().get_edgecolor() == (1, 0, 0, 1)
    assert regions[1].as_artist().get_edgecolor() == (0, 1, 0, 1)
    for i in (0, 1):
        regions[i].visual['default_style'] = None
    assert regions[0].as_artist().get_edgecolor() == (1, 0, 0, 1)
    assert regions[1].as_artist().get_edgecolor() == (0, 0, 0, 1)

    # Line2D
    assert regions[2].visual['default_style'] == 'ds9'
    assert regions[2].as_artist().get_markeredgecolor() == 'red'
    assert regions[3].as_artist().get_markeredgecolor() == '#00ff00'
    for i in (2, 3):
        regions[i].visual['default_style'] = None
    assert regions[2].as_artist().get_markeredgecolor() == 'red'
    assert regions[3].as_artist().get_markeredgecolor() == 'C0'

    # Text
    assert regions[4].visual['default_style'] == 'ds9'
    assert regions[4].as_artist().get_color() == 'red'
    assert regions[5].as_artist().get_color() == '#00ff00'
    for i in (4, 5):
        assert regions[i].as_artist().get_va() == 'center'
        assert regions[i].as_artist().get_ha() == 'center'
    for i in (4, 5):
        regions[i].visual['default_style'] = None
        assert regions[i].as_artist().get_va() == 'baseline'
        assert regions[i].as_artist().get_ha() == 'left'
    assert regions[4].as_artist().get_color() == 'red'
    assert regions[5].as_artist().get_color() == 'black'
