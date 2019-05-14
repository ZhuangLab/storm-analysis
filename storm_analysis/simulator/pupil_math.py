#!/usr/bin/env python
"""
Some math for calculating PSFs from pupil functions. All units are in 
microns.

Important note - The default for the simulator, and what is also used
in the diagnostics, is a pupil function with a pixel size of 1/2 the 
actual pixel size. This was done as it has a more realistic width. If 
you use the correct pixel size the PSF will be too narrow. This can
also be handled using OTF scaling as described in the Hanser paper,
but as this is more complicated. Also the pupil function localization
software would need to be updated to include OTF scaling, which it
currently does not.


This is based on code provided by the Huang Lab at UCSF. 

McGorty et al. "Correction of depth-dependent aberrations in 3D 
single-molecule localization and super-resolution microscopy", 
Optics Letters, 2014.


Another reference for pupil functions is:

Hanser et al. "Phase-retrieved pupil functions in wide-field fluorescence
microscopy", Journal of Microscopy, 2004.


Also, thanks to Sjoerd Stallinga for providing his MATLAB code for calculating
vectorial PSFs.

M. Siemons, C. N. Hulleman, R. Ø. Thorsen, C. S. Smith, and S. Stallinga,
"High precision wavefront control in point spread function engineering 
for single emitter localization", Optics Express, 26, pp. 8397-8416, 2018.


Hazen 05/19
"""

import math
import numpy
import scipy
import scipy.fftpack
import tifffile

import storm_analysis
import storm_analysis.pupilfn.otf_scaling_c as otfSC
import storm_analysis.pupilfn.pupil_function_c as pfFnC
import storm_analysis.simulator.pf_math_c as pfMathC


class PupilMathException(storm_analysis.SAException):
    pass


class Geometry(object):

    def __init__(self, size, pixel_size, wavelength, imm_index, NA):
        """
        size - The number of pixels in the PSF image, assumed square
               and a multiple of 2.
        pixel_size - The size of the camera pixel in um.
        wavelength - The wavelength of the flourescence in um.
        imm_index - The index of the immersion media.
        NA - The numerical aperature of the objective.
        """
        super(Geometry, self).__init__()

        if not ((size%2)==0):
            raise PupilMathException("PF size must be a multiple of 2!")

        # imm_index must be larger than the objective NA.
        if (imm_index <= NA):
            raise PupilMathException("Immersion media index must be larger than objective NA!")

        self.imm_index = float(imm_index)
        self.NA = float(NA)
        self.pixel_size = float(pixel_size)
        self.size = int(size)
        self.wavelength = float(wavelength)

        # Hanser, 2004, page 35.
        self.k_max = NA/wavelength

        dk = 1.0/(size * pixel_size)
        self.r_max = self.k_max/dk
        
        [x,y] = numpy.mgrid[ -self.size/2.0 : self.size/2.0, -self.size/2.0 : self.size/2.0]

        # Vectors to use for X/Y translation.
        self.kx = x/size
        self.ky = y/size

        kx = dk * x
        ky = dk * y
        self.k = numpy.sqrt(kx * kx + ky * ky)

        # Hanser, 2004, page 34.
        tmp = imm_index/wavelength

        # Vector to use for Z translation.
        self.kz = numpy.lib.scimath.sqrt(tmp * tmp - self.k * self.k)

        self.r = self.k/self.k_max
        
        self.kz[(self.r > 1.0)] = 0.0
        self.n_pixels = numpy.sum(self.r <= 1)
        self.norm = math.sqrt(self.r.size)

        if False:
            with tifffile.TiffWriter("kz.tif") as tf:
                tf.save(numpy.abs(self.kz).astype(numpy.float32))
                tf.save(numpy.angle(self.kz).astype(numpy.float32))

    def aberration(self, depth, smp_index):
        """
        Models the effect of a refractive index difference between the sample
        and the immersion media. See Hanser 2004, equations 4-9.
  
        depth - Point source depth in microns.
        smp_index - Refractive index of the sample media.

        Returns total aberration function (Hanser 2004, equation 8). Multiply the PF by
        this numpy array to include this aberration.
        """
        sin_theta_1 = (self.wavelength/self.imm_index)*self.k

        # Special handling of the center point where self.k = 0.0, this
        # will cause problems because then theta_1 will also be 0.0 and
        # we'll end up with 0.0/0.0 when we calculate amp_comp. So instead
        # we just use a really small number.
        #
        cp = int(sin_theta_1.shape[0]/2)
        sin_theta_1[cp,cp] = 1.0e-6
                
        sin_theta_2 = (self.imm_index/smp_index)*sin_theta_1
        
        # Limit to physical values.
        if(smp_index < self.imm_index):
            mask = (sin_theta_2 > 1.0)
        else:
            mask = (sin_theta_1 > 1.0)

        # Set values in the masked region to something that shouldn't cause 
        # problems..
        #
        sin_theta_1[mask] = 0.5
        sin_theta_2[mask] = 0.5
        
        theta_1 = numpy.arcsin(sin_theta_1)
        theta_2 = numpy.arcsin(sin_theta_2)
        
        amp_trans = (sin_theta_1 * numpy.cos(theta_2)/numpy.sin(theta_1 + theta_2))*(1.0 + (1.0/numpy.cos(theta_2 - theta_1)))
        amp_comp = (self.imm_index * numpy.tan(theta_2))/(smp_index * numpy.tan(theta_1))
        
        # This eliminates all the non-physical values.
        amp_trans[mask] = 0.0

        new_phase = numpy.pi * 2.0 * depth * (smp_index * numpy.cos(theta_2) - self.imm_index * numpy.cos(theta_1))/self.wavelength
        
        return amp_trans * amp_comp * numpy.exp(1j * new_phase)

    def aberrationOPD(self, na, ta, ne = None, te = 150.0):
        """
        Calculates an optical path difference aberration like the ones used in the
        Gibson-Lanni PSF model.

        na - Actual index of the media.

        ta - Actual media thickness (usually in microns).

        ne - Expected index of the media (default is self.imm_index).

        te - Expected media thickness (default is 150 microns).

        Returns the aberration function, multiply the PF by this numpy array to include
        this aberration.
        """
        if ne is None:
            ne = self.imm_index

        # ne must be larger than the objective NA.
        if (ne <= self.NA):
            raise PupilMathException("Expected index must be larger than objective NA!")
        
        # na must be larger than the objective NA.
        if (na <= self.NA):
            raise PupilMathException("Actual index must be larger than objective NA!")
        
        t1 = self.NA*self.r
        t1[(self.r > 1.0)] = self.NA
        t1 = t1*t1

        t2 = ne * te * numpy.sqrt(1 - t1/(ne*ne))
        t3 = na * ta * numpy.sqrt(1 - t1/(na*na))
        
        phase = numpy.pi * 2.0 * (t2 - t3) / self.wavelength
            
        pf_ab = numpy.exp(1j * phase)
        pf_ab[(self.r > 1.0)] = 0.0

        return pf_ab

    def applyNARestriction(self, pupil_fn):
        """
        pupil_fn - The pupil function to restrict the NA of.
    
        return - The NA restricted pupil function.
        """
        pupil_fn[(self.r > 1.0)] = 0.0
        return pupil_fn

    def beadScalingFactor(self, diameter):
        """
        diameter - Bead diameter in microns.

        Return a bead function to use for OTF scaling.

        Hanser, 2004, page 36.
        """
        x = numpy.pi * self.k * diameter

        # Special handling of the center point where k = 0.0, at
        # this point the correct value is 1.0.
        #
        cp = int(self.size/2)
        x[cp,cp] = 1.0e-12

        b_k = 3.0*(numpy.sin(x)/(x*x*x) - numpy.cos(x)/(x*x))
        b_k[cp,cp] = 1.0
        
        return b_k

    def changeFocus(self, pupil_fn, z_dist):
        """
        pupil_fn - The pupil function.
        z_dist - The distance to the new focal plane.
        
        return - The pupil function at the new focal plane.
        """
        return numpy.exp(1j * 2.0 * numpy.pi * self.kz * z_dist) * pupil_fn

    def createPlaneWave(self, n_photons):
        """
        n_photons - The intensity of the pupil function.
    
        return - The pupil function for a plane wave.
        """
        plane = numpy.sqrt(n_photons/self.n_pixels) * numpy.exp(1j * numpy.zeros(self.r.shape))
        if False:
            tifffile.imsave("plane.tif", numpy.angle(self.applyNARestriction(plane)).astype(numpy.float32))   
        return self.applyNARestriction(plane)

    def createFromZernike(self, n_photons, zernike_modes):
        """
        n_photons - The intensity of the pupil function
        zernike_modes - List of lists, [[magnitude (in radians), m, n], [..]]
    
        return - The pupil function for this combination of zernike modes.
        """
        if (len(zernike_modes) == 0):
            return self.createPlaneWave(n_photons)
        else:
            phases = numpy.zeros(self.r.shape)
            for zmn in zernike_modes:
                phases = pfMathC.zernikeGrid(phases, zmn[0], zmn[1], zmn[2], radius = self.r_max)
            zmnpf = numpy.sqrt(n_photons/self.n_pixels) * numpy.exp(1j * phases)
            if False:
                tifffile.imsave("zmnpf.tif", numpy.angle(self.applyNARestriction(zmnpf)).astype(numpy.float32))        
            return self.applyNARestriction(zmnpf)

    def dx(self, pupil_fn):
        """
        Returns the derivative of the pupil function in x.
        """
        return -1j * 2.0 * numpy.pi * self.kx * pupil_fn

    def gaussianScalingFactor(self, sigma):
        """
        Returns a gaussian function to use for OTF rescaling.
        """
        return numpy.exp(-1 * self.k * self.k / (2.0 * sigma * sigma))
        
    def pfToPSF(self, pf, z_vals, want_intensity = True, scaling_factor = None):
        """
        pf - A pupil function.
        z_vals - The z values (focal planes) of the desired PSF.
        want_intensity - (Optional) Return intensity, default is True.
        scaling_factor - (Optional) The OTF rescaling factor, default is None.

        return - The PSF that corresponds to pf at the requested z_vals.
        """
        if want_intensity:
            psf = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]))
            for i, z in enumerate(z_vals):
                defocused = toRealSpace(self.changeFocus(pf, z))
                if scaling_factor is not None:
                    otf = scipy.fftpack.fftshift(scipy.fftpack.fft2(intensity(defocused)))
                    otf_scaled = otf * scaling_factor
                    psf[i,:,:] = numpy.abs(scipy.fftpack.ifft2(otf_scaled))
                else:
                    psf[i,:,:] = intensity(defocused)
            return psf
        else:
            if scaling_factor is not None:
                raise PupilMathException("OTF scaling of a complex valued PSF is not supported!")
            
            psf = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                              dtype = numpy.complex_)
            for i, z in enumerate(z_vals):
                psf[i,:,:] = toRealSpace(self.changeFocus(pf, z))
            return psf

    def theoreticalOTF(self):
        """
        Returns a theoretical OTF.

        OTF = 2 * (psi - cos(psi)sin(psi)) / pi
        psi = inv_cos(lambda * wavelength / 2 * NA)

        Reference:
        https://www.microscopyu.com/microscopy-basics/modulation-transfer-function

        I'm assuming that the formula in this reference is for the FT of the PSF and 
        not the FT of the square root of the PSF.
        """
        tmp = (1.5 * self.wavelength * self.k / (2.0 * self.NA))
        tmp[(tmp > 1.0)] = 1.0
        tmp[(tmp < -1.0)] = -1.0
        
        psi = numpy.arccos(tmp)
        otf = 2.0 * (psi - numpy.cos(psi)*numpy.sin(psi)) / numpy.pi

        otf = otf/numpy.sum(otf)

        if False:
            tifffile.imsave("otf.tif", otf.astype(numpy.float32))
        
        return otf

    def translatePf(self, pupil_fn, dx, dy):
        """
        Translate the Pf using Fourier translation.
    
        pupil_fn - A pupil function.
        dx - Translation in x in pixels.
        dy - Translation in y in pixels.
    
        return - The PF translated by dx, dy.
        """
        return numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy)) * pupil_fn

    
class GeometryVectorial(Geometry):
    """
    PF model including the vectorial nature of light.
    """
    def __init__(self, size, pixel_size, wavelength, imm_index, NA):
        """
        size - The number of pixels in the PSF image, assumed square.
        pixel_size - The size of the camera pixel in um.
        wavelength - The wavelength of the flourescence in um.
        imm_index - The index of the immersion media.
        NA - The numerical aperature of the objective.
        """
        super(GeometryVectorial, self).__init__(size, pixel_size, wavelength, imm_index, NA)

        t1 = self.k*self.wavelength/self.imm_index
        t1[(t1 > 1.0)] = 1.0
        
        self.theta = numpy.arcsin(t1)
        self.theta[(self.r > 1.0)] = 0.0
        
        self.phi = numpy.arctan2(self.ky, self.kx)
        self.phi[(self.r > 1.0)] = 0.0

        c_theta = numpy.cos(self.theta)
        c_theta[(self.r > 1.0)] = 0.0

        s_theta = numpy.sin(self.theta)
        s_theta[(self.r > 1.0)] = 0.0
        
        c_phi = numpy.cos(self.phi)
        c_phi[(self.r > 1.0)] = 0.0
        
        s_phi = numpy.sin(self.phi)
        s_phi[(self.r > 1.0)] = 0.0
        
        self.px_ex = c_theta*c_phi*c_phi + s_phi*s_phi
        self.px_ey = (1.0 - c_theta)*s_phi*c_phi
        self.py_ex = (c_theta - 1.0)*s_phi*c_phi
        self.py_ey = c_theta*s_phi*s_phi + c_phi*c_phi
        self.pz_ex = s_theta*c_phi
        self.pz_ey = s_theta*s_phi
    
    def pfToPSF(self, pf, z_vals, scaling_factor = None):
        """
        This matches the super-class function of the same name, but there is no
        option of returning the (complex) real space PF. Also the PSF that is
        returned is that for a freely rotating dipole.

        pf - A pupil function.
        z_vals - The z values (focal planes) of the desired PSF.
        scaling_factor - (Optional) The OTF rescaling factor, default is None.

        return - The PSF that corresponds to pf at the requested z_vals.
        """
        rs_pf = self.pfToRS(pf, z_vals)
        psf = self.rsToPSF(rs_pf)

        if scaling_factor is not None:
            for i in range(psf.shape[0]):
                otf = scipy.fftpack.fftshift(scipy.fftpack.fft2(psf[i,:,:]))
                otf_scaled = otf * scaling_factor
                psf[i,:,:] = numpy.abs(scipy.fftpack.ifft2(otf_scaled))

        return psf
        
    def pfToRS(self, pf, z_vals):
        """
        pf - A pupil function.
        z_vals - The z values (focal planes) of the desired PSF.
 
        return - The real space counterpart of the PF at z_vals as a list.
                 [[rs_px_ex, rs_px_ey], [rs_py_ex, rs_py_ey], [rs_pz_ex, rs_pz_ey]]
        """
        rs_px_ex = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)
        rs_px_ey = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)
        
        rs_py_ex = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)
        rs_py_ey = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)

        rs_pz_ex = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)
        rs_pz_ey = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)

        for i, z in enumerate(z_vals):
            pf_at_z = self.changeFocus(pf, z)

            rs_px_ex[i,:,:] = toRealSpace(self.px_ex * pf_at_z)
            rs_px_ey[i,:,:] = toRealSpace(self.px_ey * pf_at_z)

            rs_py_ex[i,:,:] = toRealSpace(self.py_ex * pf_at_z)
            rs_py_ey[i,:,:] = toRealSpace(self.py_ey * pf_at_z)

            rs_pz_ex[i,:,:] = toRealSpace(self.pz_ex * pf_at_z)
            rs_pz_ey[i,:,:] = toRealSpace(self.pz_ey * pf_at_z)

        return [[rs_px_ex, rs_px_ey], [rs_py_ex, rs_py_ey], [rs_pz_ex, rs_pz_ey]]

    def rsToPSF(self, rs_pf, ratios = [1.0/3.0, 1.0/3.0, 1.0/3.0]):
        """
        rs_pf - List of real space PFs from the pfToRS() method.
        ratios - Ratios of dipole contribution, should sum to 1.0.
        """
        psf = numpy.zeros((rs_pf[0][0].shape[0], self.size, self.size))
        for j in range(3):
            psf += ratios[j] * intensity(rs_pf[j][0])
            psf += ratios[j] * intensity(rs_pf[j][1])

        return psf

    
class GeometryC(Geometry):
    """
    This class uses some of the C libraries in pupilfn to do the heavy lifting. It assumes
    that the OTF scaling array is symmetric in X/Y.

    Based on profiling this is 3-4x faster than the pure Python version.
    """
    def __init__(self, size, pixel_size, wavelength, imm_index, NA):
        """
        size - The number of pixels in the PSF image, assumed square.
        pixel_size - The size of the camera pixel in um.
        wavelength - The wavelength of the flourescence in um.
        imm_index - The index of the immersion media.
        NA - The numerical aperature of the objective.
        """
        super(GeometryC, self).__init__(size, pixel_size, wavelength, imm_index, NA)

        self.otf_sc = otfSC.OTFScaler(size = size)
        self.pf_c = pfFnC.PupilFunction(geometry = self)

    def __del__(self):
        self.otf_sc.cleanup()
        self.pf_c.cleanup()

    def pfToPSF(self, pf, z_vals, want_intensity = True, scaling_factor = None):
        """
        pf - A pupil function.
        z_vals - The z values (focal planes) of the desired PSF.
        want_intensity - (Optional) Return intensity, default is True.
        scaling_factor - (Optional) The OTF rescaling factor, default is None.

        return - The PSF that corresponds to pf at the requested z_vals.
        """
        self.pf_c.setPF(pf)
        
        if want_intensity:
            if scaling_factor is not None:
                self.otf_sc.setScale(scaling_factor)
                
            psf = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]))
            for i, z in enumerate(z_vals):
                self.pf_c.translateZ(z)
                temp = self.pf_c.getPSFIntensity()
                if scaling_factor is not None:
                    psf[i,:,:] = self.otf_sc.scale(temp)
                else:
                    psf[i,:,:] = temp
                
        else:
            if scaling_factor is not None:
                raise PupilMathException("OTF scaling of a complex valued PSF is not supported!")
            
            psf = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]), dtype = numpy.complex_)
            for i, z in enumerate(z_vals):
                self.pf_c.translateZ(z)
                psf[i,:,:] = self.pf_c.getPSF()
            return psf
        
        return psf
    

class GeometryCVectorial(GeometryVectorial):
    """
    The vectorial PF version of GeometryC.
    """
    def __init__(self, size, pixel_size, wavelength, imm_index, NA):
        """
        size - The number of pixels in the PSF image, assumed square.
        pixel_size - The size of the camera pixel in um.
        wavelength - The wavelength of the flourescence in um.
        imm_index - The index of the immersion media.
        NA - The numerical aperature of the objective.
        """
        super(GeometryCVectorial, self).__init__(size, pixel_size, wavelength, imm_index, NA)

        self.otf_sc = otfSC.OTFScaler(size = size)
        self.pf_c = pfFnC.PupilFunction(geometry = self)

    def __del__(self):
        self.otf_sc.cleanup()
        self.pf_c.cleanup()
    
    def pfToPSF(self, pf, z_vals, scaling_factor = None):
        """
        This matches the super-class function of the same name, but there is no
        option of returning the (complex) real space PF. Also the PSF that is
        returned is that for a freely rotating dipole.

        pf - A pupil function.
        z_vals - The z values (focal planes) of the desired PSF.
        scaling_factor - (Optional) The OTF rescaling factor, default is None.

        return - The PSF that corresponds to pf at the requested z_vals.
        """
        rs_pf = self.pfToRS(pf, z_vals)
        psf = self.rsToPSF(rs_pf)

        if scaling_factor is not None:
            self.otf_sc.setScale(scaling_factor)
            for i in range(psf.shape[0]):
                psf[i,:,:] = self.otf_sc.scale(psf[i,:,:])

        return psf
        
    def pfToRS(self, pf, z_vals):
        """
        pf - A pupil function.
        z_vals - The z values (focal planes) of the desired PSF.
 
        return - The real space counterpart of the PF at z_vals as a list.
                 [[rs_px_ex, rs_px_ey], [rs_py_ex, rs_py_ey], [rs_pz_ex, rs_pz_ey]]
        """
        rs_px_ex = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)
        rs_px_ey = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)
        
        rs_py_ex = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)
        rs_py_ey = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)

        rs_pz_ex = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)
        rs_pz_ey = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                               dtype = numpy.complex_)

        self.pf_c.setPF(pf)
        for i, z in enumerate(z_vals):
            self.pf_c.translateZ(z)
            
            rs_px_ex[i,:,:] = self.pf_c.getPXEX()
            rs_px_ey[i,:,:] = self.pf_c.getPXEY()

            rs_py_ex[i,:,:] = self.pf_c.getPYEX()
            rs_py_ey[i,:,:] = self.pf_c.getPYEY()

            rs_pz_ex[i,:,:] = self.pf_c.getPZEX()
            rs_pz_ey[i,:,:] = self.pf_c.getPZEY()

        return [[rs_px_ex, rs_px_ey], [rs_py_ex, rs_py_ey], [rs_pz_ex, rs_pz_ey]]

    def rsToPSF(self, rs_pf, ratios = [1.0/3.0, 1.0/3.0, 1.0/3.0]):
        """
        rs_pf - List of real space PFs from the pfToRS() method.
        ratios - Ratios of dipole contribution, should sum to 1.0.
        """
        psf = numpy.zeros((rs_pf[0][0].shape[0], self.size, self.size))
        for j in range(3):
            psf += ratios[j] * intensity(rs_pf[j][0])
            psf += ratios[j] * intensity(rs_pf[j][1])

        return psf
        

class GeometrySim(Geometry):
    """
    This class is used in the simulations. It divides the pixel size
    by two so that simulations look more realistic without having
    to add the overhead of OTF scaling.
    """
    def __init__(self, size, pixel_size, wavelength, imm_index, NA):
        """
        size - The number of pixels in the PSF image, assumed square.
        pixel_size - The size of the camera pixel in um.
        wavelength - The wavelength of the flourescence in um.
        imm_index - The index of the immersion media.
        NA - The numerical aperature of the objective.
        """
        super(GeometrySim, self).__init__(size, 0.5*pixel_size, wavelength, imm_index, NA)
        

def intensity(x):
    """
    x - The (numpy array) to convert to intensity.

    return - The product of x and the complex conjugate of x.
    """
    return numpy.abs(x * numpy.conj(x))


def toRealSpace(pupil_fn):
    """
    pupil_fn - A pupil function.

    return - The pupil function in real space (as opposed to fourier space).
    """
    return scipy.fftpack.ifftshift(math.sqrt(pupil_fn.size) * scipy.fftpack.ifft2(pupil_fn))


if (__name__ == "__main__"):

    import pickle
    import sys

    if (len(sys.argv) < 2):
        print("usage: <psf> <zmn.txt> <amp>")
        exit()

    pixel_size = 0.10
    #pixel_size = 0.020
    wavelength = 0.6
    refractive_index = 1.5
    numerical_aperture = 1.4
    z_range = 1.0
    z_pixel_size = 0.010

    geo = Geometry(int(20.0/pixel_size),
                   pixel_size,
                   wavelength,
                   refractive_index,
                   numerical_aperture)

    if (len(sys.argv) == 4):
        zmn = []
        amp = float(sys.argv[3])
        with open(sys.argv[2]) as fp:
            for line in fp:
                data = line.strip().split(" ")
                if (len(data) == 3):
                    zmn.append([amp * float(data[2]), int(data[0]), int(data[1])])
    else:
        zmn = [[1.3, 2, 2]]
        #zmn = [[1, 0, 4]]
        #zmn = []

    pf = geo.createFromZernike(1.0, zmn)
    z_values = numpy.arange(-z_range, z_range + 0.5 * z_pixel_size, z_pixel_size)
    psfs = geo.pfToPSF(pf, z_values)

    #xy_size = 2.0*psfs.shape[0]
    xy_size = 100
    xy_start = int(0.5 * (psfs.shape[1] - xy_size) + 1)
    xy_end = int(xy_start + xy_size)
    psfs = psfs[:,xy_start:xy_end,xy_start:xy_end]
    
    if True:
        tifffile.imsave("kz.tif", numpy.real(geo.kz).astype(numpy.float32))
            
    if False:
        tifffile.imsave("pf_abs.tif", numpy.abs(pf).astype(numpy.float32))
        tifffile.imsave("pf_angle.tif", (180.0 * numpy.angle(pf)/numpy.pi + 180).astype(numpy.float32))

    if False:
        with tifffile.TiffWriter(sys.argv[1]) as psf_tif:
            temp = (psfs/numpy.max(psfs)).astype(numpy.float32)
            psf_tif.save(temp)

    if False:
        with open("z_offset.txt", "w") as fp:
            for i in range(z_values.size):
                fp.write("1 {0:.6f}\n".format(1000.0 * z_values[i]))
        
    if False:
        psfs = (65000.0 * (psfs/numpy.max(psfs))).astype(numpy.uint16)
        psf_dict = {"pixel_size" : pixel_size,
                    "wavelength" : wavelength,
                    "refractive_index" : refractive_index,
                    "numerical_aperture" : numerical_aperture,
                    "z_range" : z_range,
                    "zmm" : zmn,
                    "psf" : psfs}
        pickle.dump(psf_dict, open(sys.argv[1], "wb"), protocol = 2)


#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#


