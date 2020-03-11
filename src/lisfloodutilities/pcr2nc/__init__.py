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
