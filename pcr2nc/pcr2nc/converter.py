"""
This module contains the class code that executes conversion by collaborating with PCRasterReader and NetCDFWriter classes.
"""
from pcr2nc.reader import PCRasterReader
from pcr2nc.writer import NetCDFWriter


def convert(config):
    input_set = config['input_set']
    reader = PCRasterReader(input_set)
    pcr_metadata = reader.get_metadata_from_set()
    writer = NetCDFWriter(config.get('output_filename') or config.get('variable'),
                          config['metadata'],
                          pcr_metadata,
                          mapstack=not reader.input_is_single_file())
    for pcr_map, time_step in reader.fileset:
        writer.add_to_stack(pcr_map, time_step)
        pcr_map.close()
    writer.finalize()
