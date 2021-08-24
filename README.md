# Overview

This script produces a ldd based on MERIT Hydro for the clonemap region at the clonemap resolution. The script is entirely "self-containing" and only requires one to specify the path of the clonemap.

The script:
1. downloads 90-m MERIT Hydro upstream area data;
1. resamples the upstream area data using a maximum filter (to retain the river network);
1. uses the resampled upstream area data as elevation; and
1. computes the ldd.

The script requires the cross-platform wget and tar utilities to download and extract the MERIT Hydro data.


