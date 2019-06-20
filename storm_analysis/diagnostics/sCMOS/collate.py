#!/usr/bin/env python
"""
Collate analysis results.

Hazen 09/17
"""
import glob

import storm_analysis.diagnostics.collate as collateResults

import storm_analysis.diagnostics.sCMOS.settings as settings


def collate():
    dirs = sorted(glob.glob("test*"))

    if(len(dirs) == 0):
        print("No test directories found.")
        return

    collateResults.collateDAO(dirs, settings)


def collateRQE():
    dirs = sorted(glob.glob("test*"))

    if(len(dirs) == 0):
        print("No test directories found.")
        return

    collateResults.collateRQE(dirs, settings)
    
    
if (__name__ == "__main__"):
    collate()
    
