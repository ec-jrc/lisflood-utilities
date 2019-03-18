#!/usr/bin/env python
from setuptools import setup, find_packages


packages_deps = ['numpy', 'pyyaml', 'netCDF4']

setup_args = dict(name='pcr2nc',
                  version='0.1',
                  description='Convert PCRaster files to netCDF4',
                  license="mit",
                  install_requires=packages_deps,
                  author="Domenico Nappo",
                  author_email="domenico.nappo@gmail.com",
                  packages=find_packages(),
                  keywords="netCDF4 PCRaster",
                  entry_points={'console_scripts': ['pcr2nc = pcr2nc_script:main_script']},
                  zip_safe=True)


setup(**setup_args)

