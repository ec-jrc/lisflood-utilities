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

from lisfloodutilities.readers import PCRasterReader
from lisfloodutilities.writers.nc import NetCDFWriter


def convert(dataset=None, output=None, metadata=None):

    reader = PCRasterReader(dataset)
    pcr_metadata = reader.get_metadata_from_set()
    metadata.update(pcr_metadata)
    writer = NetCDFWriter(
        output or metadata['variable'].get('shortname', 'pcr2nc_output'),
        is_mapstack=not reader.input_is_single_file(),
        **metadata,

    )
    for pcr_map, time_step in reader.fileset:
        writer.add_to_stack(pcr_map, time_step)
        pcr_map.close()
    writer.finalize()
