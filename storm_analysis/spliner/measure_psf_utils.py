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


class ZScaler(object):
    """
    Used in PSF measurement to convert a floating point z value into
    a z index.
    """
    def __init__(self, z_range, z_step):
        super(ZScaler, self).__init__()

        assert(z_range > 0.0), "The z range must be positive."
        assert(z_step > 0.0), "The z step must be positive."
        assert(z_range >= z_step), "The z range must be greater than or equal to the step size."
    
        # Assert that the z_step size is a multiple of the z_range.
        assert ((int(z_range*1.0e+3) % int(z_step*1.0e+3)) == 0), "The z range must be a multiple of the z step."

        self.z_mid = int(round(z_range/z_step))
        self.z_max = 2 * self.z_mid + 1
        self.z_step = z_step

    def convert(self, z):
        return int(round(z/self.z_step) + self.z_mid)

    def getMaxZ(self):
        return self.z_max
    
    def inRange(self, zi):
        return ((zi > -1) and (zi < self.z_max))


def alignPSFs(psfs, max_xy = 2, max_z = 2, max_reps = 10, verbose = True):
    """
    Align multiple PSFs in x,y,z.

    psfs - A list of PSFs, each of these has shape (nz, nxy, nxy).
    max_xy - The maximum expected alignment error xy in pixels.
    max_z - The maximum expected alignment error in z in z steps.
    max_reps - Maximum number of cycles of refinement.
    verbose - Verbose, or not.

    Returns the average PSF after alignment.
    """

    # Create working list for aligned PSFs.
    aligned_psfs = []
    for i in range(len(psfs)):
        aligned_psfs.append(psfs[i])

    starting_score = psfCorrelation(aligned_psfs)
        
    # Repeat aligning a PSF to the average of all the other PSFs.
    for i in range(max_reps):
        moving = False
        for j in range(len(psfs)):

            # Compute average of all the PSFs except the current PSF.
            sum_psf =  averagePSF(aligned_psfs, skip = j)

            # Align the current PSF to the average PSF and update
            # the list of aligned PSFs.
            #
            psf_aligner = imgCorr.Align3DProductNewtonCG(sum_psf,
                                                         xy_margin = max_xy,
                                                         z_margin = max_z)

            psf_aligner.setOtherImage(aligned_psfs[j])

            [aligned_psfs[j], q_score, disp] = psf_aligner.align()

            # Check if the PSF was translated.
            if not numpy.allclose(numpy.zeros(disp.size), disp, atol = 1.0e-3):
                moving = True
            
            if verbose:
                print(i, j, q_score, disp)

        current_score = psfCorrelation(aligned_psfs)
        
        # Print current score.
        if verbose:
            print("Quality score: {0:.6f}".format(current_score/starting_score))
            print()
            
        # Stop if the PSFs are no longer being adjusted.
        if not moving:
            break
        
        i += 1

    # Compute average of aligned PSFs.
    return [averagePSF(aligned_psfs), current_score/starting_score]


def averagePSF(psfs, skip = -1):
    """
    Compute average of a list of PSFs.
    """
    n_psfs = 0
    average_psf = numpy.zeros_like(psfs[0])
    for i in range(len(psfs)):
        if (i == skip):
            continue
        average_psf += psfs[i]
        n_psfs += 1

    return average_psf/float(n_psfs)


def extractAOI(frame, aoi_size, xf, yf, zoom = 1):
    """
    Extract AOI for PSF measurements.

    frame - An image.
    aoi_size - 1/2 the AOI size in pixels.
    xf - AOI x offset in pixels.
    yf - AOI y offset in pixels.
    zoom - Zoom factor, default is 2.0.
    """
    xi = int(xf)
    yi = int(yf)

    sx = xi - aoi_size
    ex = xi + aoi_size
    sy = yi - aoi_size
    ey = yi + aoi_size
    
    # Check that the slice is inside the image.
    assert (sx >= 0), "X position is too small ({0:d}).".format(sx)
    assert (sy >= 0), "Y position is too small ({0:d}).".format(sy)
    assert (ex <= frame.shape[0]), "X position is too large ({0:d}).".format(ex)
    assert (ey <= frame.shape[1]), "Y position is too large ({0:d}).".format(ey)

    # Slice.
    im_slice = frame[sx:ex,sy:ey]

    # Zoom and center.
    if(zoom != 1):
        im_slice_up = scipy.ndimage.interpolation.zoom(im_slice, zoom)
    else:
        im_slice_up = im_slice
        
    im_slice_up = scipy.ndimage.interpolation.shift(im_slice_up, (-zoom*(xf-xi), -zoom*(yf-yi)), mode='nearest')

    return im_slice_up

    
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

    z_sclr = ZScaler(z_range, z_step)
    z_index = numpy.zeros(z_offsets.shape[0], dtype = numpy.int) - 1
    for i in range(z_offsets.shape[0]):
        if (z_offsets[i][0] < 1.0e-6):
            continue
        zi = z_sclr.convert(z_offsets[i][1])
        if z_sclr.inRange(zi):
            z_index[i] = zi
            
            #if (z_offsets[i][1] > (-z_range - 0.5*z_step)) and (z_offsets[i][1] < (z_range + 0.5*z_step)):

    assert(numpy.max(z_index) > -0.5), "No valid frames for PSF measurement."
                
    return z_index


def meanEdge(psf_slice):
    """
    Return the mean of the boundary pixels of a PSF slice.
    """
    edge = numpy.concatenate((psf_slice[0,:],
                              psf_slice[-1,:],
                              psf_slice[:,0],
                              psf_slice[:,-1]))
    return numpy.mean(edge)


def measureSinglePSFBeads(frame_reader, z_index, aoi_size, x, y, drift_xy = None, zoom = 1):
    """
    Measures a single PSF from a PSF z stack movie that you
    might take using beads.

    frame_reader - A sa_library.analysis_io.FrameReader like object.
    z_index - Z slice in the PSF for each frame, as returned for
              example by makeZIndexArray().
    aoi_size - Size of the PSF AOI.
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

        # Extract AOI.
        im_slice_up = extractAOI(frame, aoi_size, xf, yf, zoom = zoom)

        # Update accumulators.
        zi = z_index[i]
        psf[zi,:,:] += im_slice_up
        samples[zi] += 1

    return [psf, samples]


def psfCorrelation(psfs):
    """
    Calculate the correlation score of the PSFs, this is just the
    sum of the product of all the PSFs.
    """
    product = numpy.copy(psfs[0])
    for i in range(1,len(psfs)):
        product = product * psfs[i]
    product = product/float(len(psfs))
    return numpy.sum(product)


def psfSharpness(psf):
    """
    Calculates how 'sharp' the PSF is as defined here by how large 
    the mean frequency component is. The idea is that a better average
    PSF will be less blurred out, so it will have more power in
    the larger frequencies.
    """
    psd = numpy.abs(numpy.fft.fftn(psf))**2

    k1 = numpy.abs(numpy.fft.fftfreq(psf.shape[0]))
    k2 = numpy.abs(numpy.fft.fftfreq(psf.shape[1]))
    k3 = numpy.abs(numpy.fft.fftfreq(psf.shape[2]))

    # Ignore the highest frequencies as these are mostly pixel noise.
    k1[(k1 > 0.4)] = 0
    k2[(k2 > 0.4)] = 0
    k2[(k3 > 0.4)] = 0

    [m_k1, m_k2, m_k3] = numpy.meshgrid(k1, k2, k3, indexing = 'ij')
    return numpy.mean(psd * m_k1 * m_k2 * m_k3)


def smoothPSF(psf, xy_sigma = 0.5, z_sigma = 0.5):
    """
    Apply gaussian smoothing to a PSF.
    """
    return scipy.ndimage.filters.gaussian_filter(psf,
                                                 [z_sigma, xy_sigma, xy_sigma],
                                                 mode = "nearest")

    
def sumPSF(psfs):
    """
    Compute sum of a list of PSFs.
    """
    sum_psf = numpy.zeros_like(psfs[0])
    for psf in psfs:
        sum_psf += psf

    return sum_psf
