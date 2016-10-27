#!/usr/bin/env python
# -*- coding: utf-8 -*-

import platform
import os
import sys

from setuptools import setup, find_packages
from distutils.core import Extension

import numpy

version = "1.0"
description = "Read and write image data from and to TIFF files."
long_description = ""

def get_c_extensions():

    include_dirs = [os.path.join(sys.prefix, "include")]
    library_dirs = []

    if platform.system() == 'Windows':
        include_dirs += [os.environ['LIBRARY_INC']]
        library_dirs += [os.environ['LIBRARY_LIB']]
    elif platform.system() == 'Linux':
        include_dirs += []
        library_dirs += []
    elif platform.system() == 'Darwin':
        include_dirs += []
        library_dirs += []

    extensions = [#Extension("", ["./storm_analysis/fista/fista_decon_utilities.c"], ),
                  #Extension("", ["./storm_analysis/fista/fista_fft.c"], ),
                  Extension("storm_analysis.sa_library._matched_filter", ["./storm_analysis/sa_library/matched_filter.c"],
                            libraries=library_dirs + ["fftw3"], include_dirs=include_dirs + []),
                  Extension("storm_analysis.sa_library._grid", ["./storm_analysis/sa_library/grid.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.sa_library._multi_fit", ["./storm_analysis/sa_library/multi_fit.c"],
                            libraries=library_dirs + ["lapack"], include_dirs=include_dirs + []),
                  Extension("storm_analysis.sa_library._ia_utilities", ["./storm_analysis/sa_library/ia_utilities.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  #Extension("", ["./storm_analysis/dbscan/dbscan.c"], ),
                  #Extension("", ["./storm_analysis/dbscan/kdtree.c"], ),
                  #Extension("", ["./storm_analysis/decon_storm/mlem_sparse.c"], ),
                  #Extension("", ["./storm_analysis/sCMOS/scmos_utilities.c"], ),
                  #Extension("", ["./storm_analysis/L1H/fista_lib.c"], ),
                  #Extension("", ["./storm_analysis/L1H/homotopy_storm.c"], ),
                  #Extension("", ["./storm_analysis/L1H/homotopy_sse.c"], ),
                  #Extension("", ["./storm_analysis/L1H/homotopy_general.c"], ),
                  #Extension("", ["./storm_analysis/L1H/homotopy_imagea.c"], ),
                  #Extension("", ["./storm_analysis/L1H/homotopy_common.c"], ),
                  #Extension("", ["./storm_analysis/L1H/homotopy_imagea_common.c"], ),
                  #Extension("", ["./storm_analysis/L1H/homotopy_gpu.c"], ),
                  #Extension("", ["./storm_analysis/sa_utilities/fitz.c"], ),
                  #Extension("", ["./storm_analysis/sa_utilities/tracker.c"], ),
                  #Extension("", ["./storm_analysis/sa_utilities/avemlist.c"], ),
                  #Extension("", ["./storm_analysis/sa_utilities/apply-drift-correction.c"], ),
                  #Extension("", ["./storm_analysis/frc/frc.c"], ),
                  #Extension("", ["./storm_analysis/simulator/draw_gaussians.c"], ),
                  #Extension("", ["./storm_analysis/simulator/zernike.c"], ),
                  #Extension("", ["./storm_analysis/spliner/cubic_spline.c"], ),
                  #Extension("", ["./storm_analysis/spliner/multi_fit_core.c"], ),
                  #Extension("", ["./storm_analysis/spliner/cubic_fit.c"], ),
                  #Extension("", ["./storm_analysis/rolling_ball_bgr/rolling_ball_lib.c"], ),
                  ]

    return extensions

setup(
    name='storm_analysis',
    version=version,
    description=description,
    long_description=long_description,
    author='Hazen Babcock',
    author_email='hbabcock at fas.harvard.edu',
    url='https://github.com/ZhuangLab/storm-analysis',

    zip_safe=False,
    packages=find_packages(),

    ext_modules=get_c_extensions(),
    package_data={
        #'sample': ['package_data.dat'],
        # If any package contains *.txt or *.rst files, include them:
        #'': ['*.txt', '*.rst'],
        # And include any *.msg files found in the 'hello' package, too:
        #'hello': ['*.msg'],
    },
    exclude_package_data={
        #'': ['README.txt']
    },
    include_package_data=True,

    requires=['numpy (>=1.8.2)', 'setuptools'],
    
    license="",  
    keywords='storm,microscopy',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        "Programming Language :: C",
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
    ],
)