#!/usr/bin/env python
"""
This is a Python version of the C library for testing purposes.

Hazen 10/17
"""
import numpy

class PSFFFT(object):
        
    def __init__(self, psf = None, **kwds):
        super(PSFFFT, self).__init__(**kwds)

        self.psf_shape = psf.shape
        self.z_mid = int(self.psf_shape[0]/2)

        self.kx = 2.0 * numpy.pi * numpy.fft.fftfreq(self.psf_shape[2])
        self.ky = 2.0 * numpy.pi * numpy.fft.fftfreq(self.psf_shape[1])                
        self.kz = 2.0 * numpy.pi * numpy.fft.fftfreq(self.psf_shape[0])
        
        self.psf_fft = numpy.fft.fftn(psf)
        self.ws = numpy.copy(self.psf_fft)

    def getMid(self, np_arr):
        return numpy.real(np_arr[self.z_mid,:,:])
    
    def getPSF(self):
        return self.getMid(numpy.fft.ifftn(self.ws))

    def getPSFdx(self):
        tmp = numpy.copy(self.ws)

        dk_dx = -1j * self.kx
        for i in range(self.psf_shape[0]):
            for j in range(self.psf_shape[1]):
                for k in range(self.psf_shape[2]):
                    tmp[i,j,k] = self.ws[i,j,k] * dk_dx[k]
        
        return self.getMid(numpy.fft.ifftn(tmp))

    def getPSFdy(self):
        tmp = numpy.copy(self.ws)

        dk_dy = -1j * self.ky
        for i in range(self.psf_shape[0]):
            for j in range(self.psf_shape[1]):
                for k in range(self.psf_shape[2]):
                    tmp[i,j,k] = self.ws[i,j,k] * dk_dy[j]
        
        return self.getMid(numpy.fft.ifftn(tmp))

    def getPSFdz(self):
        tmp = numpy.copy(self.ws)

        dk_dz = -1j * self.kz
        for i in range(self.psf_shape[0]):
            for j in range(self.psf_shape[1]):
                for k in range(self.psf_shape[2]):
                    tmp[i,j,k] = self.ws[i,j,k] * dk_dz[i]
        
        return -self.getMid(numpy.fft.ifftn(tmp))
     
    def translate(self, dx, dy, dz):
        """
        This not supposed to be efficient, it is supposed to match
        how the C library does this calculation.

        We use -z to match the Z convention of simulator.pupilmath
        """
        k_dx = numpy.exp(-1j * dx * self.kx)
        k_dy = numpy.exp(-1j * dy * self.ky)
        k_dz = numpy.exp(-1j * -dz * self.kz)

        for i in range(self.psf_shape[0]):
            for j in range(self.psf_shape[1]):
                for k in range(self.psf_shape[2]):
                    self.ws[i,j,k] = k_dz[i]*k_dy[j]*k_dx[k]*self.psf_fft[i,j,k]

        
