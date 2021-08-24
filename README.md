# Overview

This script produces a local drain direction (LDD) map based on MERIT Hydro for the clone map region at the clone map resolution. The script is entirely "self-containing" and only requires one to specify the path of the clonemap.

The script:
1. downloads and extracts the 90-m MERIT Hydro upstream area data;
1. resamples the upstream area data using a maximum filter (to retain the river network);
1. uses the resampled upstream area data as elevation; and
1. computes the LDD.

The script requires the cross-platform wget and tar utilities to download and extract the MERIT Hydro data.

# Instructions

First modify `config.cfg` by entering the correct folders and paths. The script can subsequently be run as follows:
```
conda create --name <env> --file requirements.txt
conda activate <env>
python create_ldd_from_merit.py
```

