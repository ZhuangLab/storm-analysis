#!/usr/bin/env python
"""
Various utility functions for PSF measurement. Basically
trying to consolidate/improve what is common between the 
several different scripts that do this.

Hazen 03/18
"""
import numpy
import scipy
import scipy.ndimage

import storm_analysis.sa_library.imagecorrelation as imgCorr


def alignPSFs(psfs, max_xy = 2, max_z = 2):
    """
    Align multiple PSFs in x,y,z.

    psfs - A list of PSFs, each of these has shape (nz, nxy, nxy).
    max_xy - The maximum expected alignment error xy in pixels.
    max_z - The maximum expected alignment error in z in z steps.
    """
    aligned_psf = numpy.copy(psfs[0])
    for i in range(1, len(psfs)):
        psf_aligner = imgCorr.Align3D(aligned_psf,
                                      xy_margin = max_xy,
                                      z_margin = max_z)
        temp = psf_aligner.align(psfs[i])
        aligned_psf += temp
                   
    return aligned_psf/float(len(psfs))

    
def makeZIndexArray(z_offsets, z_range, z_step):
    """
    Create the array that specifies which slice the image at
    a particular z offset should be added to. If the image 
    should not be added to any slice then z_index will have
    the value of -1.

    Note: The bins are centered on the z_step.

    All units are in microns.

    z_offsets - The different z offsets, an array of shape
                (N,2) as contained for example in z_offsets.txt
                file.
    z_range - The range the PSF will cover (+- z_range).
    z_step - The z step size.
    """
    assert(len(z_offsets.shape) == 2), "Z offsets must have shape (N,2)."
    assert(z_offsets.shape[1] == 2), "Z offsets must have shape (N,2)."
    assert(z_range > 0.0), "The z range must be positive."
    assert(z_step > 0.0), "The z step must be positive."
    assert(z_range >= z_step), "The z range must be greater than or equal to the step size."
    
    # Assert that the z_step size is a multiple of the z_range.
    assert ((int(z_range*1.0e+3) % int(z_step*1.0e+3)) == 0), "The z range must be a multiple of the z step."
    
    z_mid = int(z_range/z_step)
    z_index = numpy.zeros(z_offsets.shape[0], dtype = numpy.int) - 1
    for i in range(z_offsets.shape[0]):
        if (z_offsets[i][0] < 1.0e-6):
            continue
        if (z_offsets[i][1] > (-z_range - 0.5*z_step)) and (z_offsets[i][1] < (z_range + 0.5*z_step)):
            z_index[i] = int(round(z_offsets[i][1]/z_step) + z_mid)

    return z_index


def measureSinglePSFBeads(frame_reader, z_index, aoi_size, x, y, drift_xy = None, zoom = 2):
    """
    Measures a single PSF from a PSF z stack movie that you
    might take using beads.

    frame_reader - A sa_library.analysis_io.FrameReader like object.
    z_index - Z slice in the PSF for each frame, as returned for
              example by makeZIndexArray().
    xy_size - Size of the PSF AOI.
    x - Bead center position in x.
    y - Bead center position in y.
    drift_xy - An array containing x,y drift information. This should
               have a shape of (N,2). The x drift is the first entry and
               the y drift is the second entry.
    zoom - Amount to magnify the final PSF by. Must be an integer.

    Returns - [psf, samples per z section]
    """
    if drift_xy is not None:
        assert(drift_xy.shape[0] == z_index.size), "XY drift must have the same number of points a z_index."
        assert(drift_xy.shape[1] == 2), "XY drift can only have an x and a y offset for each frame."

    assert(isinstance(aoi_size, int)), "PSF AOI must be an integer."
    assert(isinstance(zoom, int)), "Zoom must be an integer."

    z_size = numpy.max(z_index) + 1
    psf = numpy.zeros((z_size, 2*aoi_size*zoom, 2*aoi_size*zoom))
    samples = numpy.zeros(z_size, dtype = numpy.int)
    for i in range(z_index.size):

        # Ignore frames with 'bad' z index.
        if(z_index[i] < 0):
            continue

        # Load the frame.
        frame = frame_reader.loadAFrame(i)

        # Figure out where to slice.
        xf = x
        yf = y

        # Apply drift correction (if specified).
        if drift_xy is not None:
            xf += drift_xy[i,0]
            yf += drift_xy[i,1]
            
        xi = int(xf)
        yi = int(yf)
        zi = z_index[i]

        # Slice.
        im_slice = frame[xi - aoi_size:xi + aoi_size,
                         yi - aoi_size:yi + aoi_size]

        # Zoom and center.
        im_slice_up = scipy.ndimage.interpolation.zoom(im_slice, zoom)
        im_slice_up = scipy.ndimage.interpolation.shift(im_slice_up, (-zoom*(xf-xi), -zoom*(yf-yi)), mode='nearest')

        # Update accumulators.
        psf[zi,:,:] += im_slice_up
        samples[zi] += 1

    return [psf, samples]
