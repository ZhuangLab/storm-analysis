#!/usr/bin/python
#
# Create A matrix for compressed sensing.
#
# Hazen 08/13
# Jeff 10/12
#

import numpy
import math
import pickle
import scipy.io
from scipy.interpolate import griddata
import os

#
# Creates x,y matrix containing coordinates for PSF centers.
#
def createCoordinates(inner_size, scale, boundary_scale, pad = [0,0]):

    # Calculate number of coordinates for keep region
    n_coords = inner_size*inner_size*scale*scale
    for [i, bound_scale] in enumerate(boundary_scale):
        n_coords += ((inner_size + 2*i)*4 + 4) * bound_scale*bound_scale

    # Calculate extra boundary pixels (at the same density as the inner scale)
    num_nonpad_pixels = 2*len(boundary_scale) + inner_size
    n_coords += 4*pad[0]*num_nonpad_pixels*pad[1] # Add pixels on the edges
    n_coords += 4*pad[1]*pad[1] #Add pixels in the corners
    num_boundary = len(boundary_scale)
    print "Generating coordinate system with ", n_coords, " elements"

    # Allocate memory
    x = numpy.zeros(n_coords)
    y = numpy.zeros(n_coords)
    area = numpy.zeros(n_coords)

    # Initialize index
    current_index = 0
    
    # Fill center, measurement region
    for i in range(inner_size*scale): #y
        for j in range(inner_size*scale): #x
            x[current_index] = float(j)/float(scale) + 0.5/float(scale)
            y[current_index] = float(i)/float(scale) + 0.5/float(scale)
            area[current_index] = 1/float(scale)*1/float(scale)
            current_index += 1

    # Fill boundary regions
    for [i, bound_scale] in enumerate(boundary_scale):
        # Fill bottom and bottom corners
        for j in range(bound_scale):
            for k in range(bound_scale*(inner_size+2*i) + bound_scale*2):
                x[current_index] = float(k)/float(bound_scale) - (i+1) + 0.5/float(bound_scale)
                y[current_index] = float(j)/float(bound_scale) - (i+1) + 0.5/float(bound_scale)
                area[current_index] = 1/float(bound_scale)*1/float(bound_scale)
                current_index += 1
                
        # Fill both sides
        for j in range(bound_scale*(inner_size+2*i)):
            for k in range(bound_scale): # left side
                x[current_index] = float(k)/float(bound_scale) - (i+1) + 0.5/float(bound_scale)
                y[current_index] = float(j)/float(bound_scale) - i + 0.5/float(bound_scale)
                area[current_index] = 1/float(bound_scale)*1/float(bound_scale)
                current_index += 1
                
            for k in range(bound_scale): # right side
                x[current_index] = float(k)/float(bound_scale) + (inner_size + i) + 0.5/float(bound_scale)
                y[current_index] = float(j)/float(bound_scale) - i + 0.5/float(bound_scale)
                area[current_index] = 1/float(bound_scale)*1/float(bound_scale)
                current_index += 1
                
        # Fill top and top corners
        for j in range(bound_scale):
            for k in range(bound_scale*(inner_size+2*i) + bound_scale*2):
                x[current_index] = float(k)/float(bound_scale) - (i+1) + 0.5/float(bound_scale)
                y[current_index] = float(j)/float(bound_scale) + (inner_size + i) + 0.5/float(bound_scale)
                area[current_index] = 1/float(bound_scale)*1/float(bound_scale)
                current_index += 1

    # Add extra padding pixels (pad[0] = number of up sampled pixels per real pixel;
    #                           pad[1] = number of rows/columns of extra pixels

    # Fill bottom and bottom corners
    for j in range(pad[1]):
        for k in range(2*pad[1] + pad[0]*num_nonpad_pixels):
            x[current_index] = float(k)/float(pad[0]) - num_boundary - pad[1]/float(pad[0]) + 0.5/float(pad[0])
            y[current_index] = float(j)/float(pad[0]) - num_boundary - pad[1]/float(pad[0]) + 0.5/float(pad[0])
            area[current_index] = 1/float(pad[0])*1/float(pad[0])
            current_index += 1

    # Fill in both sides
    for j in range(pad[0]*num_nonpad_pixels):
        for k in range(pad[1]):
            x[current_index] = float(k)/float(pad[0]) - num_boundary - pad[1]/float(pad[0]) + 0.5/float(pad[0])
            y[current_index] = float(j)/float(pad[0]) - num_boundary + 0.5/float(pad[0])
            area[current_index] = 1/float(pad[0])*1/float(pad[0])
            current_index += 1

        for k in range(pad[1]):
            x[current_index] = float(k)/float(pad[0]) + (inner_size + num_boundary) + 0.5/float(pad[0])
            y[current_index] = float(j)/float(pad[0]) - num_boundary + 0.5/float(pad[0])
            area[current_index] = 1/float(pad[0])*1/float(pad[0])
            current_index += 1

    # Fill top and top corners
    for j in range(pad[1]):
        for k in range(2*pad[1] + pad[0]*num_nonpad_pixels):
            x[current_index] = float(k)/float(pad[0]) - num_boundary - pad[1]/float(pad[0]) + 0.5/float(pad[0])
            y[current_index] = float(j)/float(pad[0]) + (inner_size + num_boundary) + 0.5/float(pad[0])
            area[current_index] = 1/float(pad[0])*1/float(pad[0])
            current_index += 1
            
    # Center coordinate system
    total_dim = inner_size
    x = x-float(total_dim)/2
    y = y-float(total_dim)/2
    return [x,y,area]

#
# Generate "image analysis" A matrix from a gaussian PSF
#
def generateAMatrix(keep_pixels, keep_scale, boundary_scale, meas_pixels, sigma=1.0,normalize=True, pad = [0,0]):

    #Create coordinate system for potential fluorophores
    [x_fluor, y_fluor, area_fluor] = createCoordinates(keep_pixels, keep_scale, boundary_scale, pad)

    #Create coordinate system for pixels
    [x_pixel, y_pixel, area_pixel] = createCoordinates(meas_pixels, 1, [0])

    # Allocate A matrix
    nrows = x_pixel.size
    ncols = x_fluor.size + 1
    A = numpy.zeros((nrows,ncols))
    # Fill A matrix
    for i in range(nrows):
        for j in range(ncols-1):
            dx = x_pixel[i] - x_fluor[j]
            dy = y_pixel[i] - y_fluor[j]
            A[i,j] = gaussianPSF(dx, dy, sigma)

    if normalize:
        norms = numpy.sum(A, axis = 0)
        for i in range(ncols-1):
            A[:,i] = A[:,i]/norms[i]
    
    # Add background term
    A[:,-1] = numpy.ones(nrows)

    return A

#
# Generate "image analysis" A matrix from a measured PSF
#
def generateAMatrixFromPSF(keep_pixels, keep_scale, boundary_scale, meas_pixels, PSF_file_name, psf_scale = 1.0, pad = [0,0]):

    #Create coordinate system for potential fluorophores
    [x_fluor, y_fluor, area_fluor] = createCoordinates(keep_pixels, keep_scale, boundary_scale, pad = pad)

    #Create coordinate system for pixels
    [x_pixel, y_pixel, area_pixel] = createCoordinates(meas_pixels, 1, [0])
    
    # Read PSF
    psf_image = numpy.load(PSF_file_name)
    [psf_x_size, psf_y_size] = psf_image.shape
    
    # Build PSF Coordinate System
    psf_coord_x = numpy.linspace(-float(psf_x_size)/2.0, float(psf_x_size)/2.0, psf_x_size)
    psf_coord_y = numpy.linspace(-float(psf_y_size)/2.0, float(psf_y_size)/2.0, psf_y_size)
    psf_coord_x = psf_coord_x/psf_scale
    psf_coord_y = psf_coord_y/psf_scale

    # Interpolate PSF
    psf_interp = scipy.interpolate.RectBivariateSpline(psf_coord_x, psf_coord_y, psf_image)

    # Build Coordinate System for Finding the Center of the PSF
    psf_upsample_x = numpy.linspace(-float(psf_x_size)/2.0, float(psf_x_size)/2.0, psf_x_size*keep_scale);
    psf_upsample_y = numpy.linspace(-float(psf_y_size)/2.0, float(psf_y_size)/2.0, psf_y_size*keep_scale);

    X_grid, Y_grid = numpy.meshgrid(psf_upsample_x,psf_upsample_y)
    
    upsample_image = psf_interp.ev(X_grid.flatten(), Y_grid.flatten())
    
    max_index = numpy.argmax(upsample_image)
    index_x, index_y = numpy.unravel_index(max_index, (len(psf_upsample_x), len(psf_upsample_y)))
    
    # Reinterpolate PSF
    psf_interp = scipy.interpolate.RectBivariateSpline(psf_coord_x-psf_upsample_x[index_x], psf_coord_y-psf_upsample_y[index_y], psf_image)
    
    # Allocate A matrix
    nrows = x_pixel.size
    ncols = x_fluor.size + 1
    A = numpy.zeros((nrows,ncols))

    # Fill A matrix
    for i in range(nrows):
        for j in range(ncols-1):
            dx = x_pixel[i] - x_fluor[j]
            dy = y_pixel[i] - y_fluor[j]
            # Interpolate
            A[i,j] = psf_interp.ev(dx, dy)

    # Normalize
    norms = numpy.sum(A, axis = 0)
    for i in range(ncols-1):
        A[:,i] = A[:,i]/norms[i]
    
    # Add background term
    A[:,-1] = numpy.ones(nrows)

    return A

#
# Generate a gaussian PSF based on dx, dy, sigma
#
def gaussianPSF(dx, dy, sigma = 1.0):
    xt = dx*dx
    yt = dy*dy
    z = 1/(2*math.pi*sigma*sigma)*numpy.exp(-(xt+yt)/(2.0*sigma*sigma))
    return z

#
# Load a A matrix
#
def loadAMatrix(file_name):
    return(pickle.load(open(file_name)))

#
# Save a A matrix
#
def saveAMatrix(file_name, a_mat, meas_pixels, keep_pixels, keep_scale):
    dict = {"a_matrix" : a_mat,
            "meas_pixels" : meas_pixels,
            "keep_pixels" : keep_pixels,
            "keep_scale" : keep_scale}

    pickle.dump(dict, open(file_name, "w"))
    print "Saved " + file_name    

#
# Save A matrix in dax format for visualization purposes.
#
def saveAsDax(file_name, A, measured_pixels):
    import sa_library.daxwriter as daxwriter

    dx = daxwriter.DaxWriter(file_name,0,0)
    ncols = A.shape[1]
    for i in range(A.shape[1]):
        x = numpy.zeros(ncols)
        x[i] = 1.0
        b = numpy.dot(A,x)
        b = b.reshape(measured_pixels,measured_pixels)
        dx.addFrame(10000.0*b)

    dx.close()

#
# Visualize A matrix coordinate system using matplotlib
#
def visualize(x, y):
    import matplotlib
    import matplotlib.pyplot as pyplot
    
    fig = pyplot.figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.grid(linestyle='-',color='0.2')
    ax.plot(x,y,'o',markeredgecolor='0.0',markerfacecolor='0.0',markersize=2.0)
    pyplot.xticks(numpy.arange(13)-6.0)
    pyplot.yticks(numpy.arange(13)-6.0)
    pyplot.show()

#
# Useful function calls
#
if __name__ == "__main__":

    import sys

    if (len(sys.argv) != 4):
        print "usage: <type (theoritical | measured)> <name> <(sigma | psf.npy)>"
        exit()

    #
    # Set A matrix parameters
    #

    # The number of pixels from the original image to be analyzed (in low resolution camera pixels)
    meas_pixels = 7
    
    # The number of final pixels kept from the reconstruction (in low resolution camera pixels)
    keep_pixels = 5

    # The number of upsampled pixels per camera pixel 
    keep_scale = 8 
    
    # The number of upsampled pixels to include in each row of additional camera pixels beyond 
    # those kept for the reconstruction (each element indicates adding one more layer of boundary pixels)
    boundary = [8]

    # Parameters for a partial pixel after the boundary. The first number determines the 
    # upsampling of the pixels and the second determines how many upsampled pixels to include 
    #pad = [0,0]
    pad = [8,4]
    
    # Generate name
    name_string = sys.argv[2] + "_a" + str(meas_pixels) + "_k" + str(keep_pixels) + "_i" + str(keep_scale)
    for i in range(len(boundary)):
        name_string = name_string + "_o" + str(boundary[i])
    name_string = name_string + "_p" + str(pad[0]) + "_" + str(pad[1]) + ".amat"

    is_good = True
    #
    # Generate A matrix using a gaussian PSF.
    #
    # Note that sigma is in units of pixels.
    #
    if (sys.argv[1] == "theoritical"): 
        A = generateAMatrix(keep_pixels, keep_scale, boundary, meas_pixels, pad = pad, sigma = float(sys.argv[3]))
        A_matrix_name = name_string
        saveAMatrix(A_matrix_name, A, meas_pixels, keep_pixels, keep_scale)
        #scipy.io.savemat(A_matrix_name + ".mat", mdict={"A":A})
        
        saveAsDax("ia_matrix.dax", A, meas_pixels)
        
        #A.tofile(A_matrix_name + ".bin")

    #
    # Generate A matrix from a measured PSF.
    #
    elif (sys.argv[1] == "measured"):
        A = generateAMatrixFromPSF(keep_pixels, keep_scale, boundary, meas_pixels, sys.argv[3], psf_scale = 2.0, pad = pad)
    
        A_matrix_name = name_string
        saveAMatrix(A_matrix_name, A, meas_pixels, keep_pixels, keep_scale)
        #scipy.io.savemat(A_matrix_name + ".mat", mdict={"A":A})

        saveAsDax("ia_matrix.dax", A, meas_pixels)

    else:
        print "First argument should be 'theoritical' or 'measured'"
        is_good = False

    if is_good:
        [x_fluor, y_fluor, area_fluor] = createCoordinates(keep_pixels, keep_scale, boundary, pad)
        visualize(x_fluor, y_fluor)

#    if 1: # Save coordinate system in matlab format
#        coord_mat_name = "coord_A_" + name_string + ".mat"

#        [x_pixel, y_pixel, area_pixel] = createCoordinates(meas_pixels, 1, [0])
#        scipy.io.savemat("temp/" + coord_mat_name, mdict={"x_fluor":x_fluor, "y_fluor":y_fluor, "x_pixel":x_pixel, "y_pixel":y_pixel})

#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
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

