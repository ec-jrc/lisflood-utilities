import os
import logging
import itertools
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from lisfloodutilities.readers import PCRasterMap


logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger()


class Comparator:
    glob_expr = None

    def __init__(self):
        self.errors = []

    def compare_dirs(self, path_a, path_b):
        path_a = Path(path_a)
        path_b = Path(path_b)
        for fa in itertools.chain(*(path_a.glob(e) for e in self.glob_expr)):
            fb = path_b.joinpath(fa.name)
            if not fb.exists():
                logger.info('skipping %s as it is not in %s', fb.name, path_b.as_posix())
                continue
            errors = self.compare_files(fa, fb)
            if errors:
                self.errors += errors
        return self.errors

    def compare_files(self, fa, fb):
        raise NotImplementedError()


class PCRComparator(Comparator):
    glob_expr = ['**/*.[0-9][0-9][0-9]', '**/*.map']

    def compare_files(self, fa, fb):
        map_a = PCRasterMap(fa)
        map_b = PCRasterMap(fb)
        err = [f'{fa} different from {fb}'] if map_a != map_b else None
        map_a.close()
        map_b.close()
        return err



class TSSComparator(Comparator):
    glob_expr = ['**/*.tss']

    def compare_files(self, fa, fb):
        with open(fa, 'rb') as fp1, open(fb, 'rb') as fp2:
            # need to skip first line in TSS as it reports settings filename
            next(fp1)
            next(fp2)
            while True:
                b1 = fp1.readline()
                b2 = fp2.readline()
                if (b1 != b2) or (not b1 and b2) or (not b2 and b1):
                    err = [f'{fa} different from {fb}']
                    return err
                if not b1:
                    return None


class NetCDFComparator(Comparator):
    glob_expr = ['**/*.nc']

    def __init__(self, mask, atol=0.05, rtol=0.1, max_perc_diff=0.2, max_perc_large_diff=0.1):
        self.mask = Dataset(mask)
        maskvar = [k for k in self.mask.variables if len(self.mask.variables[k].dimensions) == 2][0]
        self.maskarea = np.logical_not(self.mask.variables[maskvar][:, :])
        self.atol = atol
        self.rtol = rtol
        self.max_perc_large_diff = max_perc_large_diff
        self.max_perc_diff = max_perc_diff
        self.large_diff_th = self.atol * 10
        super().__init__()

    def check_vars(self, maskarea, step, vara_step, varb_step, varname, filepath):
        step = step if step is not None else '(no time)'
        filepath = os.path.basename(filepath)
        vara_step = np.ma.masked_array(vara_step, maskarea)
        varb_step = np.ma.masked_array(varb_step, maskarea)
        vara_step = np.ma.compressed(vara_step).astype('float64')
        varb_step = np.ma.compressed(varb_step).astype('float64')
        diff_values = np.ma.abs(vara_step - varb_step)
        same_values = np.ma.allclose(diff_values, np.zeros(diff_values.shape), atol=self.atol, rtol=self.rtol)
        all_ok = vara_step.size == varb_step.size and same_values
        array_ok = np.isclose(diff_values, np.zeros(diff_values.shape), atol=self.atol, rtol=self.rtol, equal_nan=True)
        different_values_size = array_ok[~array_ok].size
        if (not all_ok) and (different_values_size > 0):
            max_diff = np.ma.amax(diff_values)
            perc_wrong = different_values_size * 100 / diff_values.size
            result = np.where(diff_values >= max_diff)
            rel_diff = max_diff * 100. / np.maximum(vara_step[result], varb_step[result])
            if np.max(rel_diff) > 0.01 and (perc_wrong >= self.max_perc_diff or (perc_wrong >= self.max_perc_large_diff and max_diff > self.large_diff_th)):
                rel_diff = np.array_str(rel_diff)
                logger.error(f'Var: {varname} - STEP {step}: {perc_wrong:3.2f}% of values are different. max diff: {max_diff:3.2f} (rel diff: {rel_diff}%)')
                logger.error('\nA: %s\nB: %s', np.array_str(vara_step[result]), np.array_str(varb_step[result]))
                return f'{filepath}/{varname}/{step} - {perc_wrong:3.2f}% of different values - max diff: {max_diff:3.2f} (rel diff: {rel_diff}%)'

    def compare_files(self, file_a, file_b):
        errors = []
        nca = Dataset(file_a)
        ncb = Dataset(file_b)
        num_dims = 3 if 'time' in nca.variables else 2
        var_name = [k for k in nca.variables if len(nca.variables[k].dimensions) == num_dims][0]
        vara = nca.variables[var_name]
        varb = ncb.variables[var_name]
        if 'time' in nca.variables:
            stepsa = nca.variables['time'][:]
            stepsb = ncb.variables['time'][:]
            if len(stepsa) != len(stepsb):
                logger.error('Different number of timesteps')
                len_stepsa = len(stepsa)
                len_stepsb = len(stepsb)
                return f'Files: {file_a} vs {file_b}: different size in time axis A:{len_stepsa} B:{len_stepsb}'
            for i, step in enumerate(stepsa):
                vara_step = vara[i][:, :]
                varb_step = varb[i][:, :]
                err = self.check_vars(self.maskarea, i, vara_step, varb_step, var_name, file_a)
                if err:
                    errors.append(err)
        else:
            vara_step = vara[:, :]
            varb_step = varb[:, :]
            err = self.check_vars(self.maskarea, None, vara_step, varb_step, var_name, file_a)
            if err:
                errors.append(err)
        nca.close()
        ncb.close()
        return errors
