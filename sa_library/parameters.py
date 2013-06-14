#!/usr/bin/python
#
# Handles parsing settings xml files.
#
# Hazen 05/12
#

import copy
from xml.dom import minidom, Node

default_params = 0

class Parameters:
    # Dynamically create the class by processing the 
    # parameters xml file.
    def __init__(self, parameters_file):
        xml = minidom.parse(parameters_file)
        settings = xml.getElementsByTagName("settings").item(0)
        for node in settings.childNodes:
            if node.nodeType == Node.ELEMENT_NODE:
                # single parameter setting
                if len(node.childNodes) == 1:
                    slot = node.nodeName
                    value = node.firstChild.nodeValue
                    type = node.attributes.item(0).value
                    if type == "int":
                        setattr(self, slot, int(value))
                    elif type == "int-array":
                        text_array = value.split(",")
                        int_array = []
                        for elt in text_array:
                            int_array.append(int(elt))
                        setattr(self, slot, int_array)
                    elif type == "float":
                        setattr(self, slot, float(value))
                    elif type == "float-array":
                        text_array = value.split(",")
                        float_array = []
                        for elt in text_array:
                            float_array.append(float(elt))
                        setattr(self, slot, float_array)
                    elif type == "string-array":
                        setattr(self, slot, value.split(","))
                    else: # everything else is assumed to be a string
                        setattr(self, slot, value)
                # multiple parameter settings.
                else:
                    print "multi parameter setting unimplemented."

        self.parameters_file = parameters_file


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
