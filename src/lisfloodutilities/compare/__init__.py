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

import itertools

from nine import IS_PYTHON2

if IS_PYTHON2:
    from pathlib2 import Path
else:
    from pathlib import Path

from .. import logger


class Comparator(object):
    """

    """
    glob_expr = None

    def __init__(self, array_equal=False, for_testing=False):
        """

        """
        self.array_equal = array_equal
        self.for_testing = for_testing
        self.errors = []

    def compare_dirs(self, path_a, path_b, skip_missing=False, timestep=None):
        """
        :param path_a
        :param path_b
        :param skip_missing (bool, default False). If True, ignore files that are in path_a and not in path_b
        :param timestep (default None). If passed, comparison happens only at the defined timestep
        """
        logger.info('Comparing %s and %s %s[skip missing: %s]', path_a, path_b, '(time: %s) ' % timestep if timestep else '', skip_missing)
        path_a = Path(path_a)
        path_b = Path(path_b)
        for fa in itertools.chain(*(path_a.glob(e) for e in self.glob_expr)):
            fb = path_b.joinpath(fa.name)
            if not fb.exists():
                if skip_missing:
                    logger.info('skipping %s as it is not in %s', fb.name, path_b.as_posix())
                    continue
                else:
                    message = '{} is missing in {}'.format(fb.name, path_b.as_posix())
                    if self.for_testing:
                        assert False, message
                    else:
                        self.errors.append(message)
                        continue

            self.compare_files(fa.as_posix(), fb.as_posix(), timestep)

    def compare_files(self, fa, fb, timestep=None):
        raise NotImplementedError()
