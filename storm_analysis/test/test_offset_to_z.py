#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.spliner.offset_to_z as offsetToZ


def test_offsetToZ_1():
    off_name = storm_analysis.getPathOutputTest("test_offset_to_z.off")

    # Create storm-control format position file.
    z_data = numpy.append(numpy.array([0.0, 0.0]), numpy.arange(-1.0, 1.01, 0.1))
    z_data = numpy.append(z_data, numpy.array([0.0]))
    
    data = numpy.zeros((z_data.size,4))
    data[:,3] = z_data
    numpy.savetxt(off_name, data)

    z_data = z_data[1:]
    
    # Convert.
    off_data = offsetToZ.offsetToZ(off_name)
    
    # Check conversion.
    assert(numpy.allclose(off_data[1:-1,0], numpy.ones(z_data.size  - 2)))
    assert(numpy.allclose(off_data[:,1], z_data))

    
def test_offsetToZ_2():
    off_name = storm_analysis.getPathOutputTest("test_offset_to_z.off")

    # Create storm-control format position file.
    z_data = numpy.append(numpy.array([0.0, 0.0]), numpy.arange(-1.0, 1.01, 0.1))
    z_data = numpy.append(z_data, numpy.array([0.0]))
    
    data = numpy.zeros((z_data.size,4))
    data[:,3] = z_data
    numpy.savetxt(off_name, data)

    z_data = z_data[1:]
    
    # Convert.
    off_data = offsetToZ.offsetToZ(off_name, all_valid = True)
    
    # Check conversion.
    assert(numpy.allclose(off_data[:,0], numpy.ones(z_data.size)))
    assert(numpy.allclose(off_data[:,1], z_data))


def test_offsetToZ_3():
    off_name = storm_analysis.getPathOutputTest("test_offset_to_z.off")

    # Create storm-control format position file.
    z_data = numpy.append(numpy.array([0.0, 0.0]), numpy.arange(-1.0, 1.01, 0.1))
    z_data = numpy.append(z_data, numpy.array([0.0]))
    
    data = numpy.zeros((z_data.size,4))
    data[:,3] = z_data
    numpy.savetxt(off_name, data)

    z_data = z_data[1:]
    
    # Convert.
    dz = 0.1
    off_data = offsetToZ.offsetToZ(off_name, dz = dz, all_valid = True)
    
    # Check conversion.
    assert(numpy.allclose(off_data[:,0], numpy.ones(z_data.size)))
    assert(numpy.allclose(off_data[:,1], z_data + dz))
    
    
if (__name__ == "__main__"):
    test_offsetToZ_1()
    test_offsetToZ_2()
    test_offsetToZ_3()
    
