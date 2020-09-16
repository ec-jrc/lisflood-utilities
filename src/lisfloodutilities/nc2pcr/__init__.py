"""

Copyright 2019 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

from lisfloodutilities.readers import NetCDFMap
from lisfloodutilities.writers import PCRasterWriter


def convert(in_filename, clonemap, out_filename):

    reader = NetCDFMap(in_filename)
    writer = PCRasterWriter(clonemap)

    for varname, values in reader.data:
        print(str(values.dtype))
        writer.write(out_filename.replace('.map', '_{}.map'.format(varname)), values)

    reader.close()
    writer.close()
