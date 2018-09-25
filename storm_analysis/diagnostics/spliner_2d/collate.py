#!/usr/bin/env python
"""
Collate analysis results.

Hazen 09/18
"""
import glob

import storm_analysis.diagnostics.collate as collateResults

import storm_analysis.diagnostics.spliner_2d.settings as settings

def collate():
    dirs = sorted(glob.glob("test*"))

    if(len(dirs) == 0):
        print("No test directories found.")
        return

    collateResults.collateDAO(dirs, settings)

if (__name__ == "__main__"):
    collate()
