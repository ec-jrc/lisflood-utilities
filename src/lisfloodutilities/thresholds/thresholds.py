"""
Copyright 2019-2023 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import argparse
import os
import sys

import numpy as np
import xarray as xr
from typing import Tuple, Union, List, Literal


def lmoments(values: np.ndarray) -> np.ndarray:
    """
    Compute first 2 L-moments of dataset (on first axis)
    :param values: N-D array of values
    :return: an estimate of the first two sample L-moments
    """

    nmoments = 3

    # we need to have at least four values in order
    # to make a sample L-moments estimation
    nvalues = values.shape[0]
    if nvalues < 4:
        raise ValueError(
            "Insufficient number of values to perform sample L-moments estimation"
        )

    # sort the values into ascending order
    values = np.sort(values, axis=0)

    sums = np.zeros((nmoments, *(values.shape[1:])))

    for i in range(1, nvalues + 1):
        z = i
        term = values[i - 1]
        sums[0] = sums[0] + term
        for j in range(1, nmoments):
            z -= 1
            term = term * z
            sums[j] = sums[j] + term

    y = float(nvalues)
    z = float(nvalues)
    sums[0] = sums[0] / z
    for j in range(1, nmoments):
        y = y - 1.0
        z = z * y
        sums[j] = sums[j] / z

    k = nmoments
    p0 = -1.0
    for _ in range(2):
        ak = float(k)
        p0 = -p0
        p = p0
        temp = p * sums[0]
        for i in range(1, k):
            ai = i
            p = -p * (ak + ai - 1.0) * (ak - ai) / (ai * ai)
            temp = temp + (p * sums[i])
        sums[k - 1] = temp
        k = k - 1

    lmoments = np.zeros((2, *(values.shape[1:])))
    lmoments[0] = sums[0]
    lmoments[1] = sums[1]

    return lmoments


def gumbel_parameters_moments(
    dis: xr.DataArray,
    dim: str = 'time'
) -> xr.Dataset:
    """Calculate the scale (sigma) and location (mu) parameters of the Gumbel distribution using method of moments.

    The method of moments is a statistical technique used to estimate the parameters of a distribution by equating
    the theoretical moments of the distribution to the sample moments.

    Parameters:
    -----------
    dis: xarray.DataArray
        An xarray DataArray containing the data over which to calculate the Gumbel parameters.
        The function computes the parameters along the dimension 'dim', typically time.
    dim: string
        Dimension in 'dis' over which statistics will be computed

    Returns:
    --------
    parameters: xarray.Dataset
        An xarray Dataset containing the estimated scale (sigma) and location (mu) parameters of the Gumbel distribution.
        The dataset includes two data variables:
        - 'sigma': xarray DataArray of the estimated scale parameters
        - 'mu': xarray DataArray of the estimated location parameters
    """

    # estimate parameters
    sigma = np.sqrt(6) / np.pi * dis.std(dim=dim, ddof=1, skipna=True)
    mu = dis.mean(dim=dim, skipna=True) - np.euler_gamma * sigma

    # create dataset
    coords = {coord: values for coord, values in dis_max.coords.items() if coord != dim}
    try:
        dims = [dim for dim in list(coords) if dim not in ['lat', 'lon']]
    except:
        dims = [dim for dim in list(coords) if dim not in ['y', 'x']]
    parameters = xr.Dataset({
        'sigma': xr.DataArray(sigma, dims=dims, coords=coords),
        'mu': xr.DataArray(mu, dims=dims, coords=coords)
    })

    return parameters


def gumbel_parameters_lmoments(dis_max: xr.DataArray, dim: str = 'time') -> xr.Dataset:
    """Calculate the scale (sigma) and location (mu) parameters of the Gumbel distribution using the method of L-moments.

    Parameters:
    -----------
    dis_max: xarray.DataArray
        An xarray DataArray containing the annual maximum series over which to calculate the Gumbel parameters.
    dim: string
        Temporal dimension in 'dis_max'

    Returns:
    --------
    parameters: xarray.Dataset
        An xarray Dataset containing the estimated scale (sigma) and location (mu) parameters of the Gumbel distribution.
        The dataset includes two data variables:
        - 'sigma': xarray DataArray of the estimated scale parameters
        - 'mu': xarray DataArray of the estimated location parameters
    """

    # estimate parameters
    lambda_coef = lmoments(dis_max)
    sigma = lambda_coef[1] / np.log(2)
    mu = lambda_coef[0] - sigma * np.euler_gamma

    # create dataset
    coords = {coord: values for coord, values in dis_max.coords.items() if coord != dim}
    try:
        dims = [dim for dim in list(coords) if dim not in ['lat', 'lon']]
    except:
        dims = [dim for dim in list(coords) if dim not in ['y', 'x']]
    parameters = xr.Dataset({
        'sigma': xr.DataArray(sigma, dims=dims, coords=coords),
        'mu': xr.DataArray(mu, dims=dims, coords=coords)
    })

    return parameters


def gumbel_function(
    return_period: Union[np.ndarray, int],
    sigma: Union[np.ndarray, float],
    mu: Union[np.ndarray, float]
) -> Union[np.ndarray, float]:
    """Compute the value of a Gumbel distribution associated to certain return periods.

    This function is equivalent to scipy.stats.gumbel_r.ppf(1 - 1 / return_period, loc=mu, scale=sigma)

    Parameters:
    -----------
    return_period: integer or numpy.ndarray
        return period
    sigma: float or numpy.ndarray
        scale parameter of the Gumbel distribution
    mu: float or numpy.ndarray
        location parameter of  the Gumbel distribution

    Returns:
    --------
    float or numpy.ndarray
        Value of the Gumbel distribution associated to the return period
    """
    return mu - sigma * np.log(np.log(return_period / (return_period - 1)))


def find_main_var(ds, path):
    variable_names = [k for k in ds.variables if len(ds.variables[k].dims) == 3]
    if len(variable_names) > 1:
        raise Exception("More than one variable in dataset {}".format(path))
    elif len(variable_names) == 0:
        raise Exception("Could not find a valid variable in dataset {}".format(path))
    else:
        var_name = variable_names[0]
    return var_name


def read_discharge(in_files):
    ds = xr.open_dataset(in_files)
    var = find_main_var(ds, in_files)
    da = ds[var]
    return da


def unmask_array(mask, template, data):
    data_unmask = np.empty_like(template)
    data_unmask[...] = np.NaN
    data_unmask[mask] = data
    return data_unmask


def create_dataset(dis_max, return_periods, mask, thresholds, sigma, mu):
    print("Creating dataset")
    ds_rp = xr.Dataset(
        coords={"lat": dis_max.coords["lat"], "lon": dis_max.coords["lon"]}
    )
    for i, rp in enumerate(return_periods):
        thres = unmask_array(mask, dis_max.isel(time=0).values, thresholds[i])
        ds_rp[f"rl_{rp}"] = (["lat", "lon"], thres)

    s = unmask_array(mask, dis_max.isel(time=0).values, sigma)
    print(s.shape)
    ds_rp[f"sigma"] = (["lat", "lon"], s)
    m = unmask_array(mask, dis_max.isel(time=0).values, mu)
    print(m.shape)
    ds_rp[f"mu"] = (["lat", "lon"], m)

    print(ds_rp)

    return ds_rp


def compute_thresholds_gumbel(
    dis_max: xr.DataArray,
    return_periods: List,
    method: Literal['l-moments', 'moments'] = 'l-moments',
    dim: str = 'time'
) -> xr.Dataset:
    """Fits the parameters of the Gumbel right distribution and estimates the discharge associated with the return periods

    Parameters:
    -----------
    dis_max: xarray.DataArray
        An xarray DataArray containing the annual maximum series over which to calculate the Gumbel parameters
    return_periods: List
        Return periods for which the associated discharge will be calculated
    method: string
        Method used to estimate the parameters of the Gumbel distribution
    dim: string
        Temporal dimension in 'dis_max'

    Returns:
    --------
    xr.Dataset
        An xarray Dataset containing the estimated Gumbel distribution parameters and discharge levels for the specified return periods.
        The dataset includes:
        - 'sigma': xarray DataArray of the estimated scale parameter of the Gumbel distribution.
        - 'mu': xarray DataArray of the estimated location parameter of the Gumbel distribution.
        - 'rp_X': xarray DataArray(s) of the discharge levels corresponding to each specified return period X.
    """

    # filter points with NaN in the first value
    mask = np.isfinite(dis_max.isel({dim: 0}))
    dis_max_masked = dis_max.where(mask, drop=True)

    print("Computing Gumbel coefficients")
    if method.lower() == 'l-moments':
        parameters = gumbel_parameters_lmoments(dis_max_masked, dim=dim)
    elif method.lower() == 'moments':
        parameters = gumbel_parameters_moments(dis_max_masked, dim=dim)
    else:
        raise ValueError(f"Invalid method '{method}'. Expected 'l-moments' or 'moments'.")

    print("Computing return periods")
    thresholds = xr.Dataset({f'rp_{rp}': gumbel_function(rp, parameters.sigma, parameters.mu) for rp in return_periods})

    return xr.merge([parameters, thresholds])


def main(argv=sys.argv):
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Utility to compute the discharge return period thresholds
        using the method of L-moments.
        Thresholds computed: [1.5, 2, 5, 10, 20, 50, 100, 200, 500]
        """,
        prog=prog,
    )
    parser.add_argument(
        "-d", "--discharge", help="Input discharge files (annual maxima)"
    )
    parser.add_argument("-o", "--output", help="Output thresholds file")

    args = parser.parse_args()

    dis = read_discharge(args.discharge)
    print(dis)

    return_periods = np.array([1.5, 2, 5, 10, 20, 50, 100, 200, 500])

    thresholds = compute_thresholds_gumbel(dis, return_periods)

    thresholds.to_netcdf(args.output)


def main_script():
    sys.exit(main())


if __name__ == "__main__":
    main_script()
