#!/usr/bin/env python

import sys
import os
import glob

from distutils.core import setup
from distutils.command.install_data import install_data
from setuptools.command.install import install

import urllib.request
import re
import tarfile
import tempfile

import splotch

import subprocess

# our custom install script
# captures the install_option --stan
class InstallCommand(install):
    user_options = install.user_options

    def initialize_options(self):
        install.initialize_options(self)

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        install.run(self)

# read the long description
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)),'README.md'),encoding='utf-8') as f:
    long_description = f.read()

# read the package requirements
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)),'requirements.txt'),encoding='utf-8') as f:
    install_requires = f.read().splitlines()

setup(name='abcSplotch',
      version=splotch.__version__,
      description='Approximate Bayesian model for cell-resolution Spatial Transcriptomics data',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author=splotch.__author__,
      author_email=splotch.__email__,
      url='https://github.com/nygctech/abcSplotch',
      license=splotch.__license__,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD 3-Clause "New" or "Revised" License (BSD-3-Clause)',
          'Programming Language :: Python :: 3'],
      packages=['splotch'],
      scripts=['bin/splotch_prepare_count_files','bin/splotch_generate_input_files'],
      # could not pass binaries through the scripts argument
      install_requires=install_requires,
      cmdclass={'install': InstallCommand}
)
