#!/usr/bin/env python
#
# (Python) rolling ball background estimation.
#
# Hazen 02/16
#

import numpy
import scipy
import scipy.ndimage

#
# Rolling ball smoothing class.
#
class PyRollingBall(object):

    #
    # ball_radius is the ball radius in pixels.
    # smoothing_sigma is the sigma of the gaussian to use for image smoothing prior
    #   estimating the background using the rolling ball.
    #
    # This assumes that the ball is relatively large compared
    # to variations in the surface topology.
    #
    def __init__(self, ball_radius, smoothing_sigma):
        self.ball_radius = ball_radius
        self.smoothing_sigma = smoothing_sigma

        self.ball_size = int(round(ball_radius * 0.5))
        br = ball_radius * ball_radius
        self.ball = numpy.zeros((2*self.ball_size+1, 2*self.ball_size+1))
        for x in range(2*self.ball_size+1):
            dx = x - self.ball_size
            for y in range(2*self.ball_size+1):
                dy = y - self.ball_size
                self.ball[x,y] = br - (dx * dx + dy * dy)
        self.ball = numpy.sqrt(self.ball)
        
    def estimateBG(self, image):
        image = image.astype(numpy.float)
        sm_image = scipy.ndimage.filters.gaussian_filter(image, self.smoothing_sigma)

        ball_image = numpy.zeros(image.shape)
        
        for x in range(image.shape[0]):
            print x
            for y in range(image.shape[1]):
                min_z = sm_image[x,y]
                for bx in range(2*self.ball_size+1):
                    cx = x + bx - self.ball_size
                    if (cx < 0):
                        continue
                    if (cx >= image.shape[0]):
                        continue    
                    for by in range(2*self.ball_size+1):
                        cy = y + by - self.ball_size
                        if (cy < 0):
                            continue
                        if (cy >= image.shape[1]):
                            continue
                        cz = sm_image[cx,cy] - self.ball[bx,by]
                        if (cz < min_z):
                            min_z = cz
                ball_image[x,y] = min_z

        ball_image += self.ball_radius
        return ball_image

    def removeBG(self, image):
        return image - self.estimateBG(image)
    

