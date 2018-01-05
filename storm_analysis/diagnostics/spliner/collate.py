#!/usr/bin/env python
"""
Collate analysis results for Spliner testing.

Hazen 09/17
"""
import glob

import storm_analysis.diagnostics.collate as collate

import settings

dirs = sorted(glob.glob("test*"))

if(len(dirs) == 0):
    print("No test directories found.")
    exit()

collate.collateSpliner(dirs, settings)
