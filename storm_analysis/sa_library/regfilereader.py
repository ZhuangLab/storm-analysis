#!/usr/bin/python
#
# Class for reading ImageJ MultiStackReg Transformation files.
# This assumes that there is only one reference image and that all
# of the transforms are of the same type.
#
# Hazen 08/10
#

import math
import numpy
import os
import re


# The file reader class & helper functions

def parseline(line, scale):
    try:
        [x, y] = line.split("\t")
        return [float(x)/float(scale), float(y)/float(scale)]
    except:
        print "Error parsing:", line

class RegReader:
    def __init__(self, filename, image_scale = 4, verbose = False):
        # keep a record of the filename
        self.filename = filename
        self.verbose = verbose

        # regex for determining "source" and "target" (reference) image
        frame_re = re.compile(r'Source img: ([\d]+) Target img: ([\d]+)')
        type_re = re.compile(r'(AFFINE|RIGID_BODY)')

        # first, determine how many frames are in the file
        number_frames = 0
        reg_file = open(self.filename, "r")
        while 1:
            line = reg_file.readline()
            if not line: break
            t = type_re.match(line)
            if t:
                number_frames += 1

        self.max_frame = number_frames + 1

        self.targ_image = 0
        self.f_to_f_matrix = []
        matrix = numpy.zeros((3,3))
        matrix[0,0] = 1.0
        matrix[1,1] = 1.0
        for i in range(self.max_frame):
            self.f_to_f_matrix.append(matrix)

        # load the data
        reg_file.seek(0)
        while 1:
            line = reg_file.readline()
            if not line: break
            t = type_re.match(line)
            if t:
                # record the type of transform
                self.type = t.group(1)

                m = frame_re.match(reg_file.readline())
                if m:
                    # "source" image number
                    src_num = int(m.group(1)) - 1

                    # "target" image number
                    targ_num = int(m.group(2)) - 1
                    self.targ_image = targ_num

                    # load "source"? transform
                    source = [parseline(reg_file.readline(), image_scale),
                              parseline(reg_file.readline(), image_scale),
                              parseline(reg_file.readline(), image_scale)]

                    reg_file.readline()
                    # load "target"? transform
                    target = [parseline(reg_file.readline(), image_scale),
                              parseline(reg_file.readline(), image_scale),
                              parseline(reg_file.readline(), image_scale)]

                    # calculate coordinate transform
                    if (self.type == "RIGID_BODY"):

                        # calculate the difference in angle between the two fields
                        dsx = source[0][0] - source[1][0]
                        dsy = source[0][1] - source[1][1]
                        angle1 = math.atan2(dsy, dsx)

                        dtx = target[0][0] - target[1][0]
                        dty = target[0][1] - target[1][1]
                        angle2 = math.atan2(dty, dtx)

                        angle = angle1 - angle2

                        # calculate the displacement in x and y
                        s_fx = math.cos(angle)*source[0][0] + math.sin(angle)*source[0][1]
                        s_fy = -math.sin(angle)*source[0][0] + math.cos(angle)*source[0][1]
                        dx = target[0][0] - s_fx
                        dy = target[0][1] - s_fy

                        if self.verbose:
                            print "dx, dy, angle:"
                            print dx, dy, angle
                            print ""
                        matrix = numpy.zeros((3,3))
                        matrix[0,0] = math.cos(angle)
                        matrix[1,0] = math.sin(angle)
                        matrix[2,0] = dx
                        matrix[0,1] = -math.sin(angle)
                        matrix[1,1] = math.cos(angle)
                        matrix[2,1] = dy
                        matrix[2,2] = 1.0

                        self.f_to_f_matrix[src_num] = matrix

                    if (self.type == "AFFINE"):
                        # calculate x & y coordinate transforms
                        matrix = numpy.ones((3,3))
                        matrix[0,0] = source[0][0]
                        matrix[1,0] = source[1][0]
                        matrix[2,0] = source[2][0]

                        matrix[0,1] = source[0][1]
                        matrix[1,1] = source[1][1]
                        matrix[2,1] = source[2][1]
                        
                        inv_matrix = numpy.linalg.inv(matrix)

                        xf_vector = numpy.zeros(3)
                        xf_vector[0] = target[0][0]
                        xf_vector[1] = target[1][0]
                        xf_vector[2] = target[2][0]

                        yf_vector = numpy.zeros(3)
                        yf_vector[0] = target[0][1]
                        yf_vector[1] = target[1][1]
                        yf_vector[2] = target[2][1]

                        v_abc = numpy.dot(inv_matrix, xf_vector)
                        v_def = numpy.dot(inv_matrix, yf_vector)

                        matrix = numpy.zeros((3,3))
                        matrix[0,0] = v_abc[0]
                        matrix[1,0] = v_abc[1]
                        matrix[2,0] = v_abc[2]

                        matrix[0,1] = v_def[0]
                        matrix[1,1] = v_def[1]
                        matrix[2,1] = v_def[2]

                        matrix[2,2] = 1.0

                        self.f_to_f_matrix[src_num] = matrix


        # propogate transformations
        if 1:
            # forward:
            index = self.targ_image + 1
            if(index < self.max_frame):
                matrix = self.f_to_f_matrix[index].copy()
                index += 1
                while(index < self.max_frame):
                    if self.verbose:
                        print "forward", index
                    temp = matrix.copy()
                    #matrix = numpy.dot(temp, self.f_to_f_matrix[index])
                    matrix = numpy.dot(self.f_to_f_matrix[index], temp)
                    self.f_to_f_matrix[index] = matrix
                    index += 1

            # backward
            index = self.targ_image - 1
            if(index > 0):
                matrix = self.f_to_f_matrix[index].copy()
                index -= 1
                while(index >= 0):
                    if self.verbose:
                        print "backward", index
                    temp = matrix.copy()
                    matrix = numpy.dot(self.f_to_f_matrix[index], temp)
                    self.f_to_f_matrix[index] = matrix
                    index -= 1

        if self.verbose:
            for i in range(self.max_frame):
                print "matrix:", i, self.targ_image
                print self.f_to_f_matrix[i]
                print ""

        
            
    def getTransformMatrix(self, index):
        if self.verbose:
            print ""
            print index
            print self.f_to_f_matrix[index]
            print ""
        return self.f_to_f_matrix[index]


# Utility functions
def applyTransform(f_to_f_matrix, x, y):
    p = numpy.ones(x.shape[0])
    input = numpy.concatenate((x[:,None],
                               y[:,None],
                               p[:,None]), axis = 1)
    result = numpy.dot(input, f_to_f_matrix)
    return [result[:,0], result[:,1]]

def applyTransformNoTranslation(f_to_f_matrix, x, y):
    p = numpy.zeros(x.shape[0])
    input = numpy.concatenate((x[:,None],
                               y[:,None],
                               p[:,None]), axis = 1)
    result = numpy.dot(input, f_to_f_matrix)
    return [result[:,0], result[:,1]]


if __name__ == "__main__":
    import sys
    reg_info = RegReader(sys.argv[1], image_scale = 1, verbose = True)
    x = numpy.array([150.0])
    y = numpy.array([150.0])
    print applyTransform(reg_info.getTransformMatrix(1), x, y)


#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
