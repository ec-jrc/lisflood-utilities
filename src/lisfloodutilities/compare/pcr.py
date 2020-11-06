import numpy as np

from lisfloodutilities import logger
from lisfloodutilities.compare import Comparator
from lisfloodutilities.readers import PCRasterMap


class PCRComparator(Comparator):
    # TODO add comparison with tolerance
    glob_expr = ['**/*.[0-9][0-9][0-9]', '**/*.map']

    def compare_files(self, file_a, file_b, timestep=None):
        logger.info('Comparing %s and %s', file_a, file_b)
        map_a = PCRasterMap(file_a)
        map_b = PCRasterMap(file_b)

        if map_a != map_b:
            map_a.close()
            map_b.close()
            message = '{} different from {}'.format(file_a, file_b)
            if self.for_testing:
                assert False, message
            else:
                self.errors.append(message)
                return message
        else:
            map_a.close()
            map_b.close()
            assert True


class TSSComparator(Comparator):
    glob_expr = ['**/*.tss']

    @staticmethod
    def find_timestep(tss_file, timestep):
        found_timestep = False
        current_line = tss_file.readline()
        while not found_timestep:
            if not current_line:
                break
            try:
                current_timestep = current_line.strip().split()[0]
            except IndexError:
                current_line = tss_file.readline()
            else:
                if current_timestep.decode() == str(timestep):
                    found_timestep = True
                    break
                else:
                    current_line = tss_file.readline()
        return current_line, found_timestep

    def __init__(self, atol=0.0001, rtol=0.001,
                 array_equal=False, for_testing=True):

        super(TSSComparator, self).__init__(array_equal=array_equal, for_testing=for_testing)
        self.atol = atol
        self.rtol = rtol

    def _findline_at_timestep(self, tss_file, timestep):
        b1, found_timestep = self.find_timestep(tss_file, timestep)
        if not found_timestep:
            message = '{} not found in {}'.format(timestep, tss_file.name)
            if self.for_testing:
                assert False, message
            else:
                self.errors.append(message)
        return b1

    def compare_lines_equal(self, file_a, file_b, timestep=None):
        tss1 = open(file_a, 'rb')
        tss2 = open(file_b, 'rb')
        # skip first line in TSS as it just reports settings filename
        tss1.readline()
        tss2.readline()

        if not timestep:
            # identical TSS files
            while True:
                b1 = tss1.readline()
                b2 = tss2.readline()
                if (b1 != b2) or (not b1 and b2) or (not b2 and b1):
                    message = '{} different from {}'.format(file_a, file_b)
                    if self.for_testing:
                        assert False, message
                    else:
                        self.errors.append(message)
                        return message
                if not b1:
                    break
        else:
            # check line at a given timestep
            b1 = self._findline_at_timestep(tss1, timestep)
            b2 = self._findline_at_timestep(tss2, timestep)

            if b1 != b2:
                message = '{} different from {} for timestep {}'.format(file_a, file_b, timestep)
                if self.for_testing:
                    assert False, message
                else:
                    self.errors.append(message)
        assert True
        return self.errors

    def compare_lines_tolerance(self, file_a, file_b, timestep=None):
        tss1 = open(file_a, 'rb')
        tss2 = open(file_b, 'rb')
        # skip first lines in TSS as it just reports settings filename
        line = b''
        numline = 1
        while 'timestep' not in line.decode():
            tss1.readline()
            line = tss2.readline()
            numline += 1

        while True:
            if not timestep:
                b1 = tss1.readline().strip().split()
                b2 = tss2.readline().strip().split()
                if not b1:
                    break
                self._check_tss_line_tol(b1, b2, file_a, file_b, numline)
            else:
                b1 = self._findline_at_timestep(tss1, timestep)
                b2 = self._findline_at_timestep(tss2, timestep)

                b1 = b1.strip().split()
                b2 = b2.strip().split()
                self._check_tss_line_tol(b1, b2, file_a, file_b, numline)
                # just one line
                break

        assert True
        return self.errors

    def _check_tss_line_tol(self, b1, b2, file_a, file_b, numline):
        array1 = np.array(b1, dtype='float64')
        array2 = np.array(b2, dtype='float64')
        if not np.allclose(array1, array2, rtol=self.rtol, atol=self.atol, equal_nan=True):
            message = '{} different from {} (line: {})\n line A: {}\n line B: {}'.format(file_a, file_b, numline,
                                                                                         array1, array2)
            message += '\n {}'.format(np.abs(array1 - array2))
            if self.for_testing:
                np.testing.assert_allclose(array1, array2, rtol=self.rtol, atol=self.atol, equal_nan=True)
                assert False, message
            else:
                self.errors.append(message)

    def compare_files(self, file_a, file_b, timestep=None):
        logger.info('Comparing %s and %s', file_a, file_b)
        if self.array_equal:
            self.compare_lines_equal(file_a, file_b, timestep)
        else:
            self.compare_lines_tolerance(file_a, file_b, timestep)