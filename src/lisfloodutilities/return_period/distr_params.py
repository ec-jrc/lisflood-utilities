import argparse
import os
import sys
import time

import numpy as np
import scipy
import xarray as xr

import lmoments


def distr_params_lmom(ds,distr):

    l1, l2, l3 = lmoments.lmoments_new(ds)

    res = dict()

    if distr == "GEV":
        
        t3 = l3/l2
        z = 2/(3+t3) - np.log(2)/np.log(3)
        
        k = 7.8590 * z + 2.9554 * np.square(z)    #shape
        sigma = l2 * k / ((1 - np.exp2(-k)) * scipy.special.gamma(1+k))    #scale
        mu = l1 + sigma * (scipy.special.gamma(1 + k) - 1) / k    #location

        res['k'] = k
        
    elif distr == "Gumbel":

        sigma = l2 / np.log(2)
        mu = l1 - sigma * 0.5772

    else:
        print("Invalid distribution")
        exit()
    
    res['sigma'] = sigma
    res['mu'] = mu
    
    res['l1'] = l1
    res['l2'] = l2
    res['l3'] = l3

    return res



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



def create_dataset(tmpl, mask, params):
    
    ds = xr.Dataset(coords={"latitude": tmpl.coords["latitude"],
                            "longitude": tmpl.coords["longitude"]})
    
    for p in params.keys():
        v = unmask_array(mask, tmpl.isel(time=0).values, params[p])
        print(f"{p} shape: ",v.shape)
        ds[f"{p}"] = (["latitude", "longitude"], v)
    
    print(f'\nOutput dataset:\n{ds}\n\n')

    return ds



def compute_params(dataset,distribution):
    
    # Mask NaN
    mask = np.isfinite(dataset.isel(time=0).values)    #exclude nan
    ds_masked = dataset.values[:, mask]
 
    # Computing distribution parameters
    start = time.time()
    params = distr_params_lmom(ds_masked,distr=distribution)
    end = time.time() - start
    print(f'Computation took {end:.1f} seconds\n')
        
    # Store result
    print(f'Create dataset to store results')
    ds_params = create_dataset(dataset, mask, params)

    return ds_params



def main(argv=sys.argv):
    prog = os.path.basename(argv[0])
    parser = argparse.ArgumentParser(
        description="""
        Utility to compute the coefficient of an extreme-value distribution
        using the method of L-moments.
        """,
        prog=prog,
    )
    parser.add_argument(
        "-i", "--input", help="Input file (Acc prec YEARLY MAX)")
    parser.add_argument(
        "-o", "--output", help="Output EV distribution parameters file")
    parser.add_argument(
        "-d", "--distribution", help="Extreme value probability distribution",
        choices=["GEV","Gumbel"], default="GEV")

    args = parser.parse_args()
    distr=args.distribution
    
    #Open the input file
    pr_ym = read_discharge(args.input)
    print(f'\nInput dataset\n{pr_ym}\n\n')

    #Compute the return period
    print(f"Computing {distr} coefficients...")
    params = compute_params(pr_ym,distribution=distr)
    
    #Write to NetCDF
    params.to_netcdf(args.output)
    print('END')


def main_script():
    #sys.exit(main())
    sys.exit(main())


if __name__ == "__main__":
    main_script()
    
    
