##########################################################################
''' @package helpers.py
Some helper functions for the unit tests.
'''
import os
# Path to the scratch directory.
scratch_dir = os.environ['SCRATCHDIR']
# Filename to use for results.
result_file = scratch_dir + "/result.mtx"
result_file2 = scratch_dir + "/result2.mtx"
# Threshold value for comparing floating point values.
THRESHOLD = 1e-4
# Threshold used for checking extrapolazation.
EXTRAPTHRESHOLD = 1e-1
