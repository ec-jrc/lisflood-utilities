import argparse
import logging
import sys
import os

import numpy as np

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger()

# TODO set this from CLI (with defaults values)
atol = 0.05
rtol = 0.1
max_perc_wrong_large_diff = 0.1
max_perc_wrong = 0.2
large_diff_th = atol * 10


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)

    def add_args(self):

        self.add_argument("-a", "--dataset_a", help='path to outputh of LisFlood version A', required=True)
        self.add_argument("-b", "--dataset_b", help='path to outputh of LisFlood version B', required=True)
        self.add_argument("-m", "--maskarea", help='path to mask', required=True)


def check_vars(maskarea, step, vara_step, varb_step, varname, filepath, errors):
    same_size = vara_step.size == varb_step.size
    vara_step = np.ma.masked_array(vara_step, maskarea)
    varb_step = np.ma.masked_array(varb_step, maskarea)
    vara_step = np.ma.compressed(vara_step).astype('float64')
    varb_step = np.ma.compressed(varb_step).astype('float64')
    diff_values = np.ma.abs(vara_step - varb_step)
    same_values = np.ma.allclose(diff_values, np.zeros(diff_values.shape), atol=atol, rtol=rtol)
    all_ok = same_size and same_values
    array_ok = np.isclose(diff_values, np.zeros(diff_values.shape), atol=atol, rtol=rtol, equal_nan=True)
    different_values_size = array_ok[~array_ok].size
    if (not all_ok) and (different_values_size > 0):
        max_diff = np.ma.amax(diff_values)
        large_diff = max_diff > large_diff_th
        perc_wrong = different_values_size * 100 / diff_values.size
        result = np.where(diff_values >= max_diff)
        rel_diff = max_diff * 100. / np.maximum(vara_step[result], varb_step[result])
        if np.max(rel_diff) > 0.01 and (
                perc_wrong >= max_perc_wrong or (perc_wrong >= max_perc_wrong_large_diff and large_diff)):
            logger.error(
                'Var: {} - STEP {}: {:3.2f}% of values are different. max diff: {:3.2f} (rel diff: {}%)'.format(
                    varname, step if step else '(no time)', perc_wrong, max_diff, np.array_str(rel_diff)))
            logger.error('\nA: %s\nB: %s', np.array_str(vara_step[result]), np.array_str(varb_step[result]))
            errors.append('{}/{}/{} - {:3.2f}% of different values - max diff: {:3.2f} (rel diff: {}%)'.format(
                os.path.basename(filepath), varname, step if step else '(no time)', perc_wrong, max_diff, np.array_str(rel_diff)))
