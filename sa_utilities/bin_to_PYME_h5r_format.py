#!/usr/bin/python
#
# Converts a .bin file into the PYME h5r format. As currently written
# this loses the category (or color) information of the localizations.
#
# This requires PYME and pytables.
# PYME is available here: http://code.google.com/p/python-microscopy/
#
# Hazen 11/13
#

import numpy 
import sys
import tables

from PYME.Acquire import MetaDataHandler

import sa_library.readinsight3 as readinsight3

if (len(sys.argv)!=4):
    print "usage: <bin file> <h5r file> <pixel size (nm)>"
    exit()

# Create the h5r file.
h5ResultsFile = tables.openFile(sys.argv[2], 'w')

# Create an empty metadata section.
resultsMDH = MetaDataHandler.HDFMDHandler(h5ResultsFile)

# Create an empty events table.
class SpoolEvent(tables.IsDescription):
   EventName = tables.StringCol(32)
   Time = tables.Time64Col()
   EventDescr = tables.StringCol(256)
   
resultsEvents = h5ResultsFile.createTable(h5ResultsFile.root, 'Events', SpoolEvent,filters=tables.Filters(complevel=5, shuffle=True))

# Define the fit results data type.
fresultdtype=[('tIndex', '<i4'),
              ('fitResults', [('A', '<f4'),
                              ('x0', '<f4'),
                              ('y0', '<f4'),
                              ('z0', '<f4'),
                              ('sigma', '<f4'),
                              ('background', '<f4'),
                              ('bx', '<f4'),
                              ('by', '<f4')]),
              ('fitError', [('A', '<f4'),
                            ('x0', '<f4'),
                            ('y0', '<f4'),
                            ('z0', '<f4'),
                            ('sigma', '<f4'), 
                            ('background', '<f4'),
                            ('bx', '<f4'),
                            ('by', '<f4')]), 
              ('resultCode', '<i4'), 
              ('slicesUsed', [('x', [('start', '<i4'),('stop', '<i4'),('step', '<i4')]),
                              ('y', [('start', '<i4'),('stop', '<i4'),('step', '<i4')]),
                              ('z', [('start', '<i4'),('stop', '<i4'),('step', '<i4')])])]

# Open the bin file.
i3_reader = readinsight3.I3Reader(sys.argv[1])
i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

# Save the localizations in h5r format.
print "Saving localizations"

pix_to_nm = float(sys.argv[3])
localization_number = 0
not_created = True
while (type(i3_block) != type(False)):

   print " saving localization", localization_number

   results = numpy.zeros(len(i3_block), fresultdtype)
   for i in range(len(i3_block)):
      results[i]['tIndex'] = int(i3_block['fr'][i])
      results[i]['fitResults']['x0'] = float(i3_block['xc'][i])*pix_to_nm
      results[i]['fitResults']['y0'] = float(i3_block['yc'][i])*pix_to_nm
      results[i]['fitResults']['z0'] = float(i3_block['zc'][i])
      results[i]['fitResults']['A'] = float(i3_block['a'][i])
      results[i]['fitResults']['background'] = float(i3_block['bg'][i])
      results[i]['fitResults']['sigma'] = float(i3_block['w'][i])
      
      # Arbitrary default values..
      results[i]['fitError']['x0'] = 10.0
      results[i]['fitError']['y0'] = 10.0
      
      localization_number += 1

   if not_created:

      # Create a new table which matches the datatype of the results and populate it.
      h5ResultsFile.createTable(h5ResultsFile.root, 'FitResults', results, filters=tables.Filters(complevel=5, shuffle=True), expectedrows=500000)
      not_created = False
      
   else:
      h5ResultsFile.root.FitResults.append(results)
      
   i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

# Flush the buffers in the pytables library.
h5ResultsFile.flush()

# Close the file.
h5ResultsFile.close()

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
