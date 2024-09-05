import yaml
from pathlib import Path
from typing import Union, Dict


class Config:
    def __init__(self, config_file):
        """
        Reads the configuration from a YAML file and sets default values if not provided.

        Parameters:
        -----------
        config_file: string or pathlib.Path
            The path to the YAML configuration file.
        """
        
        # read configuration file
        with open(config_file, 'r', encoding='utf8') as ymlfile:
            config = yaml.load(ymlfile, Loader=yaml.FullLoader)
            
        # input
        self.POINTS = Path(config['input']['points'])
        self.LDD_FINE = Path(config['input']['ldd_fine'])
        self.UPSTREAM_FINE = Path(config['input']['upstream_fine'])
        self.LDD_COARSE = Path(config['input']['ldd_coarse'])
        self.UPSTREAM_COARSE = Path(config['input']['upstream_coarse'])
        
        # output
        self.SHAPE_FOLDER = Path(config.get('output_folder', './shapefiles'))
        self.SHAPE_FOLDER.mkdir(parents=True, exist_ok=True)
        
        # conditions
        self.MIN_AREA = config['conditions'].get('min_area', 10)
        self.ABS_ERROR = config['conditions'].get('abs_error', 50)
        self.PCT_ERROR = config['conditions'].get('pct_error', 1)
            
        