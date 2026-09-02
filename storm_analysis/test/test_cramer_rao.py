#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.spliner.cramer_rao as cramerRao


def test_cramer_rao_bound_3d():
    """
    CRBound3D() passed the wrong keyword to CRSplineToPSF3D(), so it could
    not be constructed at all. It is the entry point used by the scripts
    that come with the C-Spline paper, so keep it working.
    """
    spline = storm_analysis.getData("test/data/test_spliner_psf.spline")

    crb = cramerRao.CRBound3D(spline, pixel_size = 100.0)

    # z_position is in nanometers for this class, see its doc string.
    crlb = crb.calcCRBound(20, 2000, z_position = 0.0)

    # [height, x, y, z, background] variances.
    assert (crlb.size == 5)
    assert (numpy.all(numpy.isfinite(crlb)))
    assert (numpy.all(crlb[1:4] > 0.0))

    # Localization is better in x/y than in z for this PSF.
    assert (crlb[1] < crlb[3])
    assert (crlb[2] < crlb[3])


def test_cramer_rao_astigmatism():
    """
    The x and y bounds should swap as z changes sign, since the PSF this
    spline was measured from is astigmatic.

    The z values here are in nanometers, and are well inside the range the
    spline covers. Passed as microns they would land far outside it, which
    is what used to happen to the C-Spline paper scripts.
    """
    spline = storm_analysis.getData("test/data/test_spliner_psf.spline")
    crb = cramerRao.CRBound3D(spline, pixel_size = 100.0)

    below = crb.calcCRBound(20, 2000, z_position = -200.0)
    above = crb.calcCRBound(20, 2000, z_position = 200.0)

    # Tighter in y below the focus, tighter in x above it.
    assert (below[2] < below[1])
    assert (above[1] < above[2])


def test_cramer_rao_paper_z_range():
    """
    The z range the C-Spline paper scripts sweep, -400 to 400 nanometers,
    should give finite bounds that actually vary with z. Before the units
    were reconciled these came back as nan.
    """
    spline = storm_analysis.getData("test/data/test_spliner_psf.spline")
    crb = cramerRao.CRBound3D(spline, pixel_size = 100.0)

    crx = numpy.zeros(9)
    for i in range(9):
        crlb = crb.calcCRBound(100, 4000, z_position = 100.0*(i - 4))
        assert (numpy.all(numpy.isfinite(crlb)))
        crx[i] = crlb[1]

    # The x bound tightens as z increases for this PSF, it does not sit at
    # a single clamped value.
    assert (numpy.std(crx) > 0.0)
    assert (crx[-1] < crx[0])


if (__name__ == "__main__"):
    test_cramer_rao_bound_3d()
    test_cramer_rao_astigmatism()
    test_cramer_rao_paper_z_range()
