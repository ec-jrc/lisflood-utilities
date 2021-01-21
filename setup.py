"""
Copyright 2019-2020 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.


---------------------------------------------------------------------------------------------------------------------------------------
To publish a new version of this distribution (git tags and pypi package), after pushed on main branch:

python setup.py testpypi
python setup.py publish

To upload new package on PyPi Test:
twine upload --repository-url https://test.pypi.org/legacy/ dist/*

Test package install
pip install --index-url https://test.pypi.org/simple/ lisflood-utilities==0.11.7

Installation with pip:
pip install lisflood-utilities
"""

import os
import sys
import subprocess
from shutil import rmtree

from setuptools import setup, find_packages, Command

current_dir = os.path.dirname(os.path.abspath(__file__))
readme_file = os.path.join(current_dir, 'README.md')
version_file = os.path.join(current_dir, 'src/lisfloodutilities/VERSION')

with open(readme_file, 'r') as f:
    long_description = f.read()

with open(version_file, 'r') as f:
    version = f.read().strip()


def _get_gdal_version():
    try:
        p = subprocess.Popen(['gdal-config', '--version'], stdout=subprocess.PIPE)
    except FileNotFoundError:
        raise SystemError('gdal-config not found.'
                          'GDAL seems not installed. '
                          'Please, install GDAL binaries and libraries for your system '
                          'and then install the relative pip package. '
                          'If you are using a conda environment: conda install -c conda-forge gdal')
    else:
        return p.communicate()[0].splitlines()[0].decode()


gdal_version = _get_gdal_version()


class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Publish lisflood-utilities package.'
    user_options = []

    @staticmethod
    def print_console(s):
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.print_console('Removing previous builds...')
            rmtree(os.path.join(current_dir, 'dist'))
        except OSError:
            pass

        self.print_console('Building Source and Wheel (universal) distribution...')
        os.system('{0} setup.py sdist'.format(sys.executable))

        self.print_console('Uploading the package to PyPI via Twine...')
        os.system('twine upload dist/*')

        self.print_console('Pushing git tags...')
        os.system('git tag v{0}'.format(version))
        os.system('git push --tags')

        sys.exit()


class UploadCommandTest(UploadCommand):

    def run(self):
        try:
            self.print_console('Removing previous builds...')
            rmtree(os.path.join(current_dir, 'dist'))
        except OSError:
            pass

        self.print_console('Building Source and Wheel (universal) distribution...')
        os.system('{} setup.py sdist'.format(sys.executable))

        self.print_console('Uploading the package to test PyPI via Twine...')
        os.system('twine upload --repository testpypi dist/*')

        sys.exit()


setup_args = dict(
    name='lisflood-utilities',
    package_dir={'': 'src/'},
    package_data={'lisfloodutilities': ['VERSION']},
    version=version,
    packages=find_packages('src'),
    description='A set of utilities for lisflood users. '
                'pcr2nc: Convert PCRaster files to netCDF CF 1.6; '
                'nc2pcr: Convert netCDF files ot PCRaster format; '
                'cutmaps: cut netCDF files;'
                'compare: compare two set of netcdf files',
    long_description=long_description,
    long_description_content_type='text/markdown',
    setup_requires=[
            'setuptools>=41.0', 'numpy',
    ],
    install_requires=['numpy', 'pyyaml==5.3',
                      # Can create corrupted environment if using conda,
                      # Better to install GDAL manually before to install lisflood-utilities package
                      # 'GDAL=={}'.format(gdal_version),
                      'netCDF4==1.5.3', 'toolz', 'xarray==0.15.1',
                      'dask', 'pandas==0.25.1', 'nine'],
    author="Valerio Lorini, Domenico Nappo, Lorenzo Alfieri",
    author_email="valerio.lorini@ec.europa.eu,domenico.nappo@gmail.com,lorenzo.alfieri@ec.europa.eu",
    keywords=['netCDF4', 'PCRaster', 'mapstack', 'lisflood', 'efas', 'glofas', 'ecmwf', 'copernicus'],
    license='EUPL 1.2',
    url='https://github.com/ec-jrc/lisflood-utilities',
    scripts=['bin/pcr2nc', 'bin/cutmaps', 'bin/compare', 'bin/nc2pcr', ],
    zip_safe=True,
    classifiers=[
          # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'Intended Audience :: Education',
          'Intended Audience :: Other Audience',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)',
          'Operating System :: Unix',
          'Operating System :: POSIX',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: MacOS :: MacOS X',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Physics',
    ],
    # setup.py publish to pypi.
    cmdclass={
        'publish': UploadCommand,
        'upload': UploadCommand,
        'testpypi': UploadCommandTest,
    },
)

setup(**setup_args)
