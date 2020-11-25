import datetime
import os
import uuid
from typing import Iterable

import numpy as np
from netCDF4 import Dataset, date2index
from nine import IS_PYTHON2

from ..readers import PCRasterMap

if IS_PYTHON2:
    from pathlib2 import Path
else:
    from pathlib import Path

from .. import logger
from ..writers import NetCDFWriter
from . import Comparator


class NetCDFComparator(Comparator):
    glob_expr = ['**/*.nc']

    def write_diff_files(self, filepath, varname, step, vara_step, varb_step, diff_values, lats, lons, time):
        self.diff_timesteps.append(time[step])
        filename, _ = os.path.splitext(os.path.basename(filepath))
        filepath_a = self.diff_folder.joinpath(filename + '_a.nc')
        filepath_b = self.diff_folder.joinpath(filename + '_b.nc')
        filepath_diff = self.diff_folder.joinpath(filename + '_diff.nc')
        metadata = {'variable': {'shortname': varname, 'compression': 9, 'least_significant_digit': 5},
                    'dtype': vara_step.dtype,
                    'rows': lats.size, 'cols': lons.size,
                    'lats': lats, 'lons': lons,
                    }
        if time:
            units = time.units if hasattr(time, 'units') else ''
            metadata.update({'time': {
                'units': units,
                }}
            )
        writer_a = NetCDFWriter(filepath_a.as_posix(), is_mapstack=step is not None, **metadata)
        writer_a.add_to_stack(vara_step, time_step=step)
        writer_a.finalize(timesteps=self.diff_timesteps)
        writer_b = NetCDFWriter(filepath_b.as_posix(), is_mapstack=step is not None, **metadata)
        writer_b.add_to_stack(varb_step, time_step=step)
        writer_b.finalize(timesteps=self.diff_timesteps)
        writer_diff = NetCDFWriter(filepath_diff.as_posix(), is_mapstack=step is not None, **metadata)
        writer_diff.add_to_stack(diff_values, time_step=step)
        writer_diff.finalize(timesteps=self.diff_timesteps)

    def __init__(self, mask=None, atol=0.0001, rtol=0.001,
                 max_perc_diff=0.2, max_perc_large_diff=0.1,
                 array_equal=False, for_testing=True,
                 save_diff_files=False):
        """

        """

        super(NetCDFComparator, self).__init__(array_equal=array_equal, for_testing=for_testing)
        self.maskarea = None
        self.masklats = None
        self.masklons = None
        self.y_min, self.y_max = None, None
        self.x_min, self.x_max = None, None
        if mask is not None:
            if not Path(mask).exists():
                raise ValueError('Mask map file not existing %s' % mask)
            # mask supposed to be a path to a file .nc or .map
            if mask.endswith('.nc'):
                mask = Dataset(mask)
                maskvar = [k for k in mask.variables if len(mask.variables[k].dimensions) == 2][0]
                self.maskarea = np.logical_not(mask.variables[maskvar][:, :])
                self.masklons = mask.variables['x'][:] if 'x' in mask.variables else mask.variables['lon'][:]
                self.masklats = mask.variables['y'][:] if 'y' in mask.variables else mask.variables['lat'][:]
            elif mask.endswith('.map'):
                maskmap = PCRasterMap(mask)
                self.maskarea = np.logical_not(maskmap.data)
                self.masklats = maskmap.lats
                self.masklons = maskmap.lons
            self.y_min, self.y_max = float(np.min(self.masklats)), float(np.max(self.masklats))
            self.x_min, self.x_max = float(np.min(self.masklons)), float(np.max(self.masklons))

        # tolerance thresholds
        self.atol = atol
        self.rtol = rtol
        self.max_perc_large_diff = max_perc_large_diff
        self.max_perc_diff = max_perc_diff
        self.large_diff_th = self.atol * 10
        self.save_diff_files = save_diff_files

        # diff files option
        self.diff_folder = None
        self.diff_timesteps = []
        if save_diff_files and not self.for_testing:
            self.diff_folder = self.create_diff_folder()
            logger.info('Diff files will be saved. You can find it in %s', self.diff_folder)
        elif save_diff_files and self.for_testing:
            logger.warning('You set both for_testing=True and save_diff_files=True. '
                           'These options are not compatible and diff files will not be saved')

    def compare_files(self, file_a, file_b, timestep=None):
        if timestep and isinstance(timestep, datetime.datetime):
            timestep = [timestep]
        if timestep and not isinstance(timestep, (datetime.datetime, Iterable)):
            raise ValueError('timestep must be of type datetime.datetime or a range of dates, but type {} was found'.format(str(type(timestep))))

        logger.info('Comparing %s and %s %s', file_a, file_b, '(from %s to %s)' % (min(timestep), max(timestep)) if timestep else '')

        with Dataset(file_a) as nca, Dataset(file_b) as ncb:
            num_dims = 3 if 'time' in nca.variables else 2
            var_name = [k for k in nca.variables if len(nca.variables[k].dimensions) == num_dims][0]
            vara = nca.variables[var_name]
            lons_a = nca.variables['x'][:] if 'x' in nca.variables else nca.variables['lon'][:]
            lats_a = nca.variables['y'][:] if 'y' in nca.variables else nca.variables['lat'][:]
            lons_b = ncb.variables['x'][:] if 'x' in ncb.variables else ncb.variables['lon'][:]
            lats_b = ncb.variables['y'][:] if 'y' in ncb.variables else ncb.variables['lat'][:]
            try:
                varb = ncb.variables[var_name]
            except KeyError:
                var_nameb = [k for k in ncb.variables if len(ncb.variables[k].dimensions) == num_dims][0]
                message = 'Files: {} vs {} have different variables names A:{} B:{}'.format(file_a, file_b, var_name, var_nameb)
                if self.for_testing:
                    assert False, message
                else:
                    self.errors.append(message)
                    return message
            if 'time' in nca.variables:
                if not timestep:
                    stepsa = nca.variables['time'][:]
                    stepsb = ncb.variables['time'][:]
                    if len(stepsa) != len(stepsb):
                        message = 'Files: {} vs {}: different number of steps A:{} B:{}'.format(file_a, file_b, len(stepsa), len(stepsb))
                        if self.for_testing:
                            assert False, message
                        else:
                            self.errors.append(message)
                            return message
                    for step, _ in enumerate(stepsa):
                        values_a = vara[step][:, :]
                        values_b = varb[step][:, :]
                        if not self.array_equal:
                            self.compare_arrays(values_a, values_b, var_name, step, file_a, lats_a, lons_a, lats_b, lons_b, time=nca.variables['time'])
                        else:
                            self.compare_arrays_equal(values_a, values_b, var_name, step, file_a, lats_a, lons_a, lats_b, lons_b, time=nca.variables['time'])
                else:
                    # check arrays at a given timestep
                    if not isinstance(timestep, Iterable):
                        timestep = [timestep]
                    for ts in timestep:
                        logger.info('Checking at timestep %s', ts)
                        ia = date2index(ts, nca.variables['time'], nca.variables['time'].calendar)
                        ib = date2index(ts, ncb.variables['time'], ncb.variables['time'].calendar)
                        values_a = vara[ia][:, :]
                        values_b = varb[ib][:, :]
                        if not self.array_equal:
                            self.compare_arrays(values_a, values_b, var_name, ia, file_a, lats_a, lons_a, lats_b, lons_b, time=nca.variables['time'])
                        else:
                            self.compare_arrays_equal(values_a, values_b, var_name, ia, file_a, lats_a, lons_a, lats_b, lons_b, time=nca.variables['time'])
            else:
                values_a = vara[:, :]
                values_b = varb[:, :]
                if not self.array_equal:
                    self.compare_arrays(values_a, values_b, var_name, filepath=file_a, lats_a=lats_a, lons_a=lons_a, lats_b=lats_b, lons_b=lons_b)
                else:
                    self.compare_arrays_equal(values_a, values_b, var_name, filepath=file_a, lats_a=lats_a, lons_a=lons_a, lats_b=lats_b, lons_b=lons_b)

    def compare_arrays(self, vara_step, varb_step, varname=None, step=None, filepath=None,
                       lats_a=None, lons_a=None, lats_b=None, lons_b=None, time=None):

        if self.maskarea is not None:
            vara_step, varb_step = self.cut_variables(lats_a, lats_b, lons_a, lons_b, vara_step, varb_step)
            diff_values = np.ma.abs(vara_step - varb_step)
            diff_values_no_nan = diff_values[~np.isnan(diff_values)]
            same_values = np.ma.allclose(diff_values_no_nan, np.zeros(diff_values_no_nan.shape), atol=self.atol, rtol=self.rtol)
        else:
            diff_values = np.abs(vara_step - varb_step)
            diff_values_no_nan = diff_values[~np.isnan(diff_values)]
            same_values = np.allclose(diff_values_no_nan, np.zeros(diff_values_no_nan.shape), atol=self.atol, rtol=self.rtol)

        all_ok = vara_step.size == varb_step.size and same_values
        array_ok = np.isclose(diff_values_no_nan, np.zeros(diff_values_no_nan.shape), atol=self.atol, rtol=self.rtol, equal_nan=True)
        different_values_size = array_ok[~array_ok].size

        if (not all_ok) and (different_values_size > 0):
            max_diff = np.ma.amax(diff_values_no_nan)  # returns a scalar
            perc_wrong = different_values_size * 100 / vara_step.size
            if perc_wrong >= self.max_perc_diff or (perc_wrong >= self.max_perc_large_diff and max_diff > self.large_diff_th):
                step = step if step is not None else '(no time)'
                filepath = os.path.basename(filepath) if filepath else '<mem>'
                varname = varname or '<unknown var>'
                message = '{}/{}@{} - {:3.2f}% of different values - max diff: {:3.6f}'.format(filepath, varname, step, perc_wrong, max_diff)
                logger.error(message)
                if self.for_testing:
                    assert False, message
                else:
                    self.errors.append(message)
                    if self.save_diff_files:
                        self.write_diff_files(filepath, varname, step, vara_step, varb_step, diff_values, lats_a, lons_a, time)

    def compare_arrays_equal(self, vara_step, varb_step, varname=None, step=None, filepath=None, lats_a=None, lons_a=None,
                             lats_b=None, lons_b=None, time=None):
        filepath = os.path.basename(filepath) if filepath else '<mem>'
        message = '{}/{}@{} is different'.format(filepath, varname, step)
        if self.maskarea is not None:
            vara_step, varb_step = self.cut_variables(lats_a, lats_b, lons_a, lons_b, vara_step, varb_step)

        if self.for_testing:
            np.testing.assert_array_equal(vara_step, varb_step, message)
        else:
            if not np.array_equal(vara_step, varb_step):
                self.errors.append(message)

                if self.save_diff_files:
                    if self.maskarea is not None:
                        diff_values = np.ma.abs(vara_step - varb_step)
                    else:
                        diff_values = np.abs(vara_step, varb_step)
                    self.write_diff_files(filepath, varname, step, vara_step, varb_step, diff_values, lats_a, lons_a, time)

    def cut_variables(self, lats_a, lats_b, lons_a, lons_b, vara_step, varb_step):
        # find indices
        if lats_a is not None:
            yas = np.where((lats_a >= self.y_min) & (lats_a <= self.y_max))[0][:, np.newaxis]
            xas = np.where((lons_a >= self.x_min) & (lons_a <= self.x_max))[0][np.newaxis, :]
            vara_step = vara_step[yas, xas]
        if lats_b is not None:
            ybs = np.where((lats_b >= self.y_min) & (lats_b <= self.y_max))[0][:, np.newaxis]
            xbs = np.where((lons_b >= self.x_min) & (lons_b <= self.x_max))[0][np.newaxis, :]
            varb_step = varb_step[ybs, xbs]
        vara_step = np.ma.masked_array(vara_step, self.maskarea).astype('float64')
        varb_step = np.ma.masked_array(varb_step, self.maskarea).astype('float64')
        return vara_step, varb_step

    @staticmethod
    def create_diff_folder():
        """
        :return: path of diff files
        :rtype: Path
        """
        diff_folder = Path(os.getcwd()).joinpath('./diffs/')
        diff_folder = diff_folder.joinpath(str(uuid.uuid4()).split('-')[-1])
        if not diff_folder.exists():
            diff_folder.mkdir(parents=True)
        return diff_folder
