#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Custom commands following:
#  https://seasonofcode.com/posts/how-to-add-custom-build-steps-and-commands-to-setuppy.html
#

import distutils.cmd
import platform
import os
import setuptools.command.build_py
import subprocess
import sys

from setuptools import setup, find_packages
from distutils.core import Extension

import numpy


version = "2.2"
description = "STORM movie analysis code."
long_description = ""

class SConsCommand(distutils.cmd.Command):
    """
    Custom command to run scons (http://scons.org/) to build the C libraries.
    """
    description = 'run scons to build C libraries'
    user_options = [
        ('scons-exe=', None, 'location of the scons executable'),
        ('compiler=', None, 'which C compiler to use, e.g. "mingw", ..')
    ]

    def initialize_options(self):
        self.scons_exe = ''
        self.compiler = ''

    def finalize_options(self):
        if self.scons_exe:
            assert os.path.exists(self.scons_exe), ("scon executable " + self.scons_exe + " not found")

    def run(self):
        if self.scons_exe:
            command = [self.scons_exe]
        else:
            command = ['scons']
        if self.compiler:
            command.extend(['-Q', 'compiler=' + self.compiler])

        self.announce('Running command: ' + str(command))
        try:
            subprocess.check_call(command)
        except OSError:
            print("Failed to build C libraries, is scons installed?")


setup(
    name='storm_analysis',
    version=version,
    description=description,
    long_description=long_description,
    author='Hazen Babcock',
    author_email='hbabcock at fas.harvard.edu',
    url='https://github.com/ZhuangLab/storm-analysis',

    cmdclass={
        'build_c' : SConsCommand,
    },

    zip_safe=False,
    packages=find_packages(),

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
        'Programming Language :: Python :: 3.6',
    ],
)
