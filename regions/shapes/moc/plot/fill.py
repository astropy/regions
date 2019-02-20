import numpy as np

from astropy_healpix import HEALPix
from astropy.coordinates import ICRS

from astropy.wcs.utils import skycoord_to_pixel

from .utils import build_plotting_moc
from . import culling_backfacing_cells

from astropy_healpix import level_to_nside

def compute_healpix_vertices(depth, ipix, wcs, origin):
    from matplotlib.path import Path

    path_vertices = np.array([])
    codes = np.array([])

    depth = int(depth)

    step = 1
    if depth < 3:
        step = 2

    nside = level_to_nside(depth)
    hp = HEALPix(order="nested", nside=nside, frame=ICRS())
    ipix_boundaries = hp.boundaries_skycoord(ipix, step=step)
    # Projection on the given WCS
    xp, yp = skycoord_to_pixel(ipix_boundaries, wcs=wcs)

    # Add offset origin given by the user (default is origin=(0, 0))
    xp += origin[0]
    yp += origin[1]

    c1 = np.vstack((xp[:, 0], yp[:, 0])).T
    c2 = np.vstack((xp[:, 1], yp[:, 1])).T
    c3 = np.vstack((xp[:, 2], yp[:, 2])).T
    c4 = np.vstack((xp[:, 3], yp[:, 3])).T

    if depth < 3:
        c5 = np.vstack((xp[:, 4], yp[:, 4])).T
        c6 = np.vstack((xp[:, 5], yp[:, 5])).T
        c7 = np.vstack((xp[:, 6], yp[:, 6])).T
        c8 = np.vstack((xp[:, 7], yp[:, 7])).T

        cells = np.hstack((c1, c2, c3, c4, c5, c6, c7, c8, np.zeros((c1.shape[0], 2))))

        path_vertices = cells.reshape((9*c1.shape[0], 2))
        single_code = np.array([Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])
    else:
        cells = np.hstack((c1, c2, c3, c4, np.zeros((c1.shape[0], 2))))

        path_vertices = cells.reshape((5*c1.shape[0], 2))
        single_code = np.array([Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])

    codes = np.tile(single_code, c1.shape[0])

    return path_vertices, codes

def compute_the_patches(moc, wcs, origin=(0, 0)):
    depth_ipix_d = moc.serialize(format="json")
    depth_ipix_clean_d = culling_backfacing_cells.from_moc(depth_ipix_d=depth_ipix_d, wcs=wcs)

    patches = []
    for depth, ipix in depth_ipix_clean_d.items():
        patch = compute_healpix_vertices(depth=depth,
                    ipix=ipix,
                    wcs=wcs,
                    origin=origin)
        patches.append(patch)

    return depth_ipix_clean_d.keys(), patches

def build_mpl_pathpatch(patches, wcs, **kw_mpl_pathpatch):
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch

    first_patch = patches[0]
    vertices_first_patch, codes_first_patch = first_patch
    path_vertices = np.array(vertices_first_patch)
    path_codes = np.array(codes_first_patch)

    for vertices, codes in patches[1:]:
        path_vertices = np.vstack((path_vertices, vertices))
        path_codes = np.hstack((path_codes, codes))

    path = Path(path_vertices, path_codes)
    pathpatch_mpl = PathPatch(path, **kw_mpl_pathpatch)

    return pathpatch_mpl

def fill(moc, wcs, origin, **kw_mpl_pathpatch):
    # Simplify the MOC for plotting purposes:
    # 1. Degrade the MOC if the FOV is enough big so that we cannot see the smallest HEALPix cells.
    # 2. For small FOVs, plot the MOC & POLYGONAL_MOC_FROM_FOV.
    moc_to_plot = build_plotting_moc(moc=moc, wcs=wcs)

    # If the FOV contains no cells, then moc_to_plot (i.e. the intersection between the moc
    # and the MOC created from the FOV polygon) will be empty.
    # If it is the case, we exit the method without doing anything.
    if moc_to_plot.empty():
        return None

    _, patches = compute_the_patches(moc=moc_to_plot, wcs=wcs, origin=origin)
    pathpatch_mpl = build_mpl_pathpatch(patches=patches, wcs=wcs, **kw_mpl_pathpatch)
    return pathpatch_mpl