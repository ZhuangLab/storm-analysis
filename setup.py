#!/usr/bin/env python
# -*- coding: utf-8 -*-

import platform
import os
import sys

from setuptools import setup, find_packages
from distutils.core import Extension

import numpy

version = "1.0"
description = "Zhuang lab STORM movie analysis code."
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

    extensions = [Extension("storm_analysis.fista._fista_decon_utilities", ["./storm_analysis/fista/fista_decon_utilities.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.fista._fista_fft", ["./storm_analysis/fista/fista_fft.c"],
                            libraries=library_dirs + ["fftw3"], include_dirs=include_dirs), 

                  Extension("storm_analysis.sa_library._matched_filter", ["./storm_analysis/sa_library/matched_filter.c"],
                            libraries=library_dirs + ["fftw3"], include_dirs=include_dirs),
                  Extension("storm_analysis.sa_library._grid", ["./storm_analysis/sa_library/grid.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.sa_library._dao_fit", ["./storm_analysis/sa_library/dao_fit.c",
                                                                   "./storm_analysis/sa_library/multi_fit.c"],
                            libraries=library_dirs + ["lapack"], include_dirs=include_dirs),
                  Extension("storm_analysis.sa_library._ia_utilities", ["./storm_analysis/sa_library/ia_utilities.c"],
                            libraries=library_dirs, include_dirs=include_dirs),

                  Extension("storm_analysis.dbscan._dbscan", ["./storm_analysis/dbscan/kdtree.c",
                                                              "./storm_analysis/dbscan/dbscan.c"],
                            libraries=library_dirs, include_dirs=include_dirs),

                  Extension("storm_analysis.decon_storm._mlem_sparse", ["./storm_analysis/decon_storm/mlem_sparse.c"],
                            libraries=library_dirs, include_dirs=include_dirs),

                  Extension("storm_analysis.sCMOS._scmos_utilities", ["./storm_analysis/sCMOS/scmos_utilities.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
  
                  Extension("storm_analysis.frc._frc", ["./storm_analysis/frc/frc.c"],
                            libraries=library_dirs, include_dirs=include_dirs),

                  Extension("storm_analysis.simulator._draw_gaussians", ["./storm_analysis/simulator/draw_gaussians.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.simulator._zernike", ["./storm_analysis/simulator/zernike.c"],
                            libraries=library_dirs, include_dirs=include_dirs),

                  Extension("storm_analysis.spliner._cubic_spline", ["./storm_analysis/spliner/cubic_spline.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.spliner._cubic_fit", ["./storm_analysis/spliner/cubic_fit.c",
                                                                  "./storm_analysis/sa_library/multi_fit.c",
                                                                  "./storm_analysis/spliner/cubic_spline.c"],
                            libraries=library_dirs + ["lapack"], include_dirs=include_dirs),

                  Extension("storm_analysis.rolling_ball_bgr._rolling_ball_lib", ["./storm_analysis/rolling_ball_bgr/rolling_ball_lib.c"],
                            libraries=library_dirs, include_dirs=include_dirs),

                  Extension("storm_analysis.sa_utilities._fitz", ["./storm_analysis/sa_utilities/fitz.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.sa_utilities._tracker", ["./storm_analysis/sa_utilities/tracker.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.sa_utilities._avemlist", ["./storm_analysis/sa_utilities/avemlist.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.sa_utilities._apply-drift-correction", ["./storm_analysis/sa_utilities/apply-drift-correction.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  
                  Extension("storm_analysis.L1H._homotopy_storm", ["./storm_analysis/L1H/homotopy_storm.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.L1H._homotopy_sse", ["./storm_analysis/L1H/homotopy_sse.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.L1H._homotopy_general", ["./storm_analysis/L1H/homotopy_general.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.L1H._homotopy_imagea", ["./storm_analysis/L1H/homotopy_imagea.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.L1H._homotopy_common", ["./storm_analysis/L1H/homotopy_common.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.L1H._homotopy_imagea_common", ["./storm_analysis/L1H/homotopy_imagea_common.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.L1H._homotopy_gpu", ["./storm_analysis/L1H/homotopy_gpu.c"],
                            libraries=library_dirs, include_dirs=include_dirs),
                  Extension("storm_analysis.L1H._homotopy_ia_storm", ["./storm_analysis/L1H/homotopy_imagea.c",
                                                                      "./storm_analysis/L1H/homotopy_storm.c",
                                                                      "./storm_analysis/L1H/homotopy_imagea_common.c",
                                                                      "./storm_analysis/L1H/homotopy_common.c"],
                            libraries=library_dirs + ["lapack"], include_dirs=include_dirs),
                  Extension("storm_analysis.L1H._homotopy_ia_sse", ["./storm_analysis/L1H/homotopy_imagea.c",
                                                                    "./storm_analysis/L1H/homotopy_sse.c",
                                                                    "./storm_analysis/L1H/homotopy_imagea_common.c",
                                                                    "./storm_analysis/L1H/homotopy_common.c"],
                            libraries=library_dirs + ["lapack"], include_dirs=include_dirs),
                  ]
    if platform.system() == 'Windows':
        extensions += [Extension("storm_analysis.L1H._fista_lib", ["./storm_analysis/L1H/fista_lib.c"],
                                 libraries=library_dirs, include_dirs=include_dirs),
 
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
    package_data={},
    exclude_package_data={},
    include_package_data=True,

    requires=[],

    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    
    license="",  
    keywords='storm,microscopy',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: Mixed',
        "Programming Language :: C",
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
