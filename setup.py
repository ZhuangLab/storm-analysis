#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.core import Extension

import numpy

version = "1.0"
description = "Read and write image data from and to TIFF files."
long_description = ""

def get_c_extensions():
    extensions = [#Extension("", ["./storm_analysis/fista/fista_decon_utilities.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/fista/fista_fft.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/sa_library/matched_filter.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/sa_library/grid.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  Extension("", ["./storm_analysis/sa_library/multi_fit.c"],
                            include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/sa_library/ia_utilities.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/dbscan/dbscan.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/dbscan/kdtree.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/decon_storm/mlem_sparse.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/sCMOS/scmos_utilities.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/L1H/fista_lib.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/L1H/homotopy_storm.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/L1H/homotopy_sse.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/L1H/homotopy_general.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/L1H/homotopy_imagea.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/L1H/homotopy_common.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/L1H/homotopy_imagea_common.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/L1H/homotopy_gpu.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/sa_utilities/fitz.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/sa_utilities/tracker.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/sa_utilities/avemlist.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/sa_utilities/apply-drift-correction.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/frc/frc.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/simulator/draw_gaussians.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/simulator/zernike.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/spliner/cubic_spline.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/spliner/multi_fit_core.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/spliner/cubic_fit.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
                  #Extension("", ["./storm_analysis/rolling_ball_bgr/rolling_ball_lib.c"],
                  #          include_dirs=[], library_dirs=[], extra_link_args=[], extra_compile_args=[]),
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