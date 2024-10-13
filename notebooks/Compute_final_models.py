#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 22:27:26 2024

@author: ciaran

Code to compute some example final models e.g. median and mean

Run separately for DGRF, IGRF and SV

"""

# Import the packages and bespoke functions from src folder
import os
import sys
import glob

import numpy as np
sys.path.append('..')


from src import shc_utils as shau

field_type = 'sv'   # 'main' or 'sv'
candidate = 'SV'    # 'DGRF', 'IGRF' or 'SV'
year = '2025-2030'         # For labelling

_DIR = os.path.abspath('../data/coefficients/' + candidate + '/*.cof')


# %%
"""

Check files are correctly formatted

"""
print('List of ' + candidate + ' Candidates \n ----------- ')
files = sorted(glob.glob(_DIR))
print(*files, sep='\n')

pass_or_fail = shau.check_cof_file_format(files, field_type)

if not(pass_or_fail).all():
    raise Exception('A candidate model format is incorrect. Fix it before \
                    continuing')

# Load in the coefficients - define the type (main or sv) and the leading
# acronym on the filename
[coeffs, institute_name, degree] = \
         shau.load_coeffs_to_numpy(sorted(files), field_type, candidate)
num_candidates = coeffs.shape[1]


# %% models - mean and median

xgrf_mean = np.mean(coeffs,axis=1)

# Note median is the average of the two middle values if N is even
xgrf_median = np.median(coeffs,axis=1)



# %% Write out files 
decimal_places = 2

filename = os.path.abspath('../data/coefficients/' + candidate + '/' 
                           + candidate + '_Mean.cof')
header = [candidate + ' ' + year + ' Mean']

shau.write_coefficients_to_cof(filename, xgrf_mean, [], 
                              field_type, ' '.join(header), decimal_places)


filename = os.path.abspath('../data/coefficients/' + candidate + '/' 
                           + candidate +  '_Median.cof')
header = [candidate + ' ' + year + ' Median']

shau.write_coefficients_to_cof(filename, xgrf_median, [], 
                              field_type, ' '.join(header), decimal_places)
