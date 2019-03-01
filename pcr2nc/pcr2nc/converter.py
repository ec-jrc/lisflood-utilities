"""
This module contains the class code that executes conversion by collaborating with PCRasterReader and NetCDFWriter classes.
"""
from pcr2nc.reader import PCRasterReader
from pcr2nc.writer import NetCDFWriter


class Converter:
    """
    A class representing a Converter object that is initialized with a config dictionary.
    The config dictionary contains input arguments and metadata values.
    """
    def __init__(self, config):
        self.config = config
        input_set = config['input_set']
        self.reader = PCRasterReader(input_set)
        pcr_metadata = self.reader.get_metadata_from_set()
        self.writer = NetCDFWriter(config.get('output_filename') or config.get('variable'),
                                   config['metadata'],
                                   pcr_metadata, mapstack=not self.reader.input_is_single_file())

    def convert(self):
        """
        Execute the conversion and write the NetCDF4 file.
        """
        for pcr_map, time_step in self.reader.fileset:
            self.writer.add_to_stack(pcr_map, time_step)
            pcr_map.close()
        self.writer.finalize()
