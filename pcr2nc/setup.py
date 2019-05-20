"""
Copyright 2019 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.


---------------------------------------------------------------------------------------------------------------------------------------
python setup.py sdist

To upload new package on PyPi Test:
twine upload --repository-url https://test.pypi.org/legacy/ dist/*

To upload new package on PyPi:
twine upload dist/*

Test package install
pip install --index-url https://test.pypi.org/simple/ pcr2nc==0.1

In prod:
pip install pcr2nc
"""

import os
from setuptools import setup, find_packages


packages_deps = ['numpy>=1.15', 'pyyaml>4.2b1', 'netCDF4>=1.3.1']
current_dir = os.path.dirname(os.path.abspath(__file__))
readme_file = os.path.join(current_dir, 'README.md')
version_file = os.path.join(current_dir, 'VERSION')

with open(readme_file, 'r') as f:
    long_description = f.read()

with open(version_file, 'r') as f:
    version = f.read()

setup_args = dict(name='pcr2nc',
                  version=version,
                  packages=find_packages(),
                  description='Convert PCRaster files to netCDF4',
                  long_description=long_description,
                  long_description_content_type='text/markdown',
                  install_requires=packages_deps,
                  author="Domenico Nappo",
                  author_email="domenico.nappo@gmail.com",
                  keywords=['netCDF4', 'PCRaster', 'mapstack'],
                  license='EUPL 1.2',
                  url='https://github.com/ec-jrc/lisflood-model',
                  entry_points={'console_scripts': ['pcr2nc = pcr2nc_script:main_script']},
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
                  )


setup(**setup_args)

