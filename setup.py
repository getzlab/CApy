from setuptools import setup
import re
import os
import numpy as np
import sys
from distutils.core import Extension

ver_info = sys.version_info
if ver_info < (3,5,4):
    raise RuntimeError("CApy requires at least python 3.5.4")

setup(
    name = 'capy',
    version = "0.1",
    packages = [
        'capy',
    ],
    description = 'Cancer Analysis Tools in Python',
    url = 'https://github.com/broadinstitute/capy',
    author = 'Julian Hess',
    author_email = 'jhess@broadinstitute.org',
    install_requires = [
        'pandas==1.4.2',
        'numpy>=1.16.4',
        'scipy>=1.5.4',
        'statsmodels>=0.12.2',
        'pyfaidx>=0.5.5.2',
        'tqdm>=4.32.1',
        'firecloud-dalmatian>=0.0.17',
        'matplotlib>=3.1.1',
        'adjustText>=0.7.3',
        'pyarrow>=0.17.0'
    ],
	ext_modules = [Extension(
	  'fastmmap',
	  sources = ['capy/fastmmap.c'],
	  include_dirs = [np.get_include()]
	)]
)
