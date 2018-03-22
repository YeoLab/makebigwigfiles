#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='makebigwigfiles',
    version='0.0.1',
    url='github.com/YeoLab/makebigwigfiles',
    license='',
    author='gpratt, adomissy, byee4',
    author_email='bay001@ucsd.edu',
    description='Makes strand-specific bigwig files from a BAM file',
    packages=['makebigwigfiles'],
    package_dir={
        'makebigwigfiles': 'makebigwigfiles',
    },
    entry_points = {
        'console_scripts': [
            'makebigwigfiles = makebigwigfiles.make_bigwig_files:main',
        ]
    }
)
