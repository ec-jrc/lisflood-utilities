import logging
import yaml
from pathlib import Path
from typing import Optional

# set logger
logger = logging.getLogger(__name__)

class Config:
    """
    Manages the application's configuration by reading a YAML file
    and setting default values.
    """
    
    def __init__(self, config_file: Path):
        """
        Reads the configuration from a YAML file and sets default values if not provided.

        Parameters:
        -----------
        config_file: string or pathlib.Path
            The path to the YAML configuration file.
        """
        
        # extract the working directory
        config_file = Path(config_file).resolve()
        self.base_path = config_file.parent

        # read configuration file
        with open(config_file, 'r', encoding='utf8') as ymlfile:
            config = yaml.load(ymlfile, Loader=yaml.FullLoader)
            
        # input file paths
        inputs = config['input']
        self.points: Optional[Path] = self._absolute_path(inputs.get('points', None))
        self.points_fine: Optional[Path] = self._absolute_path(inputs.get('points_fine', None))
        self.basins_fine: Optional[Path] = self._absolute_path(inputs.get('basins_fine', None))
        self.flwdir_fine: Optional[Path] = self._absolute_path(inputs.get('flwdir_fine', None))
        self.flwdir_coarse: Optional[Path] = self._absolute_path(inputs.get('flwdir_coarse', None))

        # tasks to be done
        if self.points is None:
            self.run_fine = False
            if (self.points_fine is None) or (self.basins_fine is None):
                raise ValueError("If 'points' is not provided, both 'points_fine' and 'basins_fine' need to be provided.")
            else:
                self.run_coarse = True
        else:
            self.run_fine = True
            if self.flwdir_coarse is None:
                self.run_coarse = False
            else:
                self.run_coarse = True

        # resolutions
        self.fine_resolution: Optional[str] = None
        self.coarse_resolution: Optional[str] = None
        
        # output folder
        self.output_folder: Optional[Path] = self._absolute_path(config.get('output_folder', 'output'))
        if self.output_folder is not None:
            self.output_folder.mkdir(parents=True, exist_ok=True)
        else:
            raise ValueError("Output folder path is not defined in the configuration.")
        
        # conditions
        conditions = config['conditions']
        self.min_area = conditions.get('min_area', 10)
        self.abs_error = conditions.get('abs_error', 50)
        self.pct_error = conditions.get('pct_error', 1)

    def _absolute_path(self, path_str: Optional[str]) -> Optional[Path]:
        """
        Helper function to join the "base_path" with relative paths from the YAML configuration file.
        """
        if path_str is None:
            return None
        path = Path(path_str)
        return path if path.is_absolute() else (self.base_path / path).resolve()