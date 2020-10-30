"""

Copyright 2019-2020 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import os
import sysconfig
import sys
import logging

IS_PYTHON2 = sys.version_info[0] == 2
pkg_dir = os.path.join(sysconfig.get_paths()['purelib'], 'lisfloodutilities')
version_file = os.path.join(pkg_dir, 'VERSION') if os.path.exists(pkg_dir) else './src/lisfloodutilities/VERSION'
version = open(version_file).read().strip() if os.path.exists(version_file) else '0.0.0.dev'

logging.basicConfig(format='[%(asctime)s] - %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
logger = logging.getLogger()

__version__ = []
for numb in version.split('.'):
    try:
        __version__.append(int(numb))
    except ValueError:
        # avoid to add to __version__ suffixes like 'post1'
        pass
__version__ = tuple(__version__)
