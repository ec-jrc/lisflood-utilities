"""
Copyright 2019-2020 European Union
Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");
You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:
https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt
Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""

import argparse
import sys
import numpy as np
import pcraster as pcr
import os
import time
import logging

logging.basicConfig(format='[%(asctime)s] - %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
logger = logging.getLogger()


def main(cliargs):
    parser = ParserHelpOnError(description='Define Water Regions consistent with calibration points: {}')
    parser.add_args()
    args = parser.parse_args(cliargs)

    calib_points = args.calib_points
    countries_id = args.countries_id
    ldd = args.ldd
    waterregions_initial = args.waterregions_initial
    output_wr = args.output_wr

    [waterregion, subcat1] = define_waterregions(calib_points, countries_id, ldd, waterregions_initial,output_wr)
    logger.info('\nUsing %s and %s to define the water regions\n ', calib_points, ldd)
    
    # check the consistency between the water regions and the calibration catchments
    cal_catch = pcr.pcr_as_numpy(subcat1)
    cal_catch[cal_catch < 1] = -9999
    cal_catch[np.isnan(cal_catch) == 1] = -9999
    
    wr = pcr.pcr_as_numpy(waterregion)
    wr_id = np.unique(wr)
       
    id_error_wr = []
    cal_catch_error_wr = []
    output_message = []
  
    for a in wr_id:
           extract_wr = np.where(wr == a,cal_catch,-9999)
           num_cal_catch = np.unique(extract_wr)
           num_cal_catch_check = np.extract(num_cal_catch != -9999,num_cal_catch)
           if len(num_cal_catch_check) > 1:
             a = int(a)
             id_error_wr.append(a)
             num_cal_catch_write = num_cal_catch[num_cal_catch != -9999].astype(int)
             cal_catch_error_wr.append(num_cal_catch_write)
             output_message = 'ERROR: The  water regions WR are included in more than one calibration catchment \n WR=' + str(id_error_wr) + '\n Calibration catchments= ' + str(cal_catch_error_wr)
    if id_error_wr == []:
             output_message = 'OK! Each water region is completely included inside one calibration catchment.'
    
    print(output_message)
    os.unlink('./points.map')


class ParserHelpOnError(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(1)
        
    def add_args(self):

        self.add_argument("-p", "--calib_points",
                          help='list of calibration points: lon or x, lat or y, point id',
                          required=True)
        self.add_argument("-l", "--ldd",
                          help='LDD map pcraster format',
                          required=True)
        self.add_argument("-C", "--countries_id",
                          help='map of Countries ID pcraster format',
                          required=True)
        self.add_argument("-w", "--waterregions_initial",
                          help='initial map of water regions pcraster format',
                          required=True)
        self.add_argument("-o", "--output_wr",
                          help='output map of water regions pcraster format',
                          required=True)

def define_waterregions(calib_points=None, countries_id=None, ldd=None, waterregions_initial=None, output_wr=None):

    
    #1. The calibration points are converted into a map
    command_string = 'col2map -N ' + calib_points + ' points.map --large --clone ' + countries_id
    os.system(command_string)

    #2. Cacthment map1 derived from calibration points
    ldd1 = pcr.readmap(ldd)
    points = pcr.readmap('points.map')
    subcat1 = pcr.subcatchment(ldd1,points)
   
    #3. Making map with all valid ldd cells
    land = pcr.scalar(subcat1)
    land = pcr.nominal(land >= 0.99999)
  
    #4. Taking only the catchments from calibration points
    subcat2 = pcr.ifthen(subcat1 != 0,subcat1) 

    #5. Splitting riverbasins if part of different countries, so if cross-border catchments
    scalar_a = pcr.scalar(subcat2)
    scalar_b = pcr.scalar(countries_id)
    subcat3 = pcr.nominal(scalar_a*scalar_b)

    #6. Covering the missing areas with the old waterregions
    subcat4 = pcr.cover(subcat3,waterregions_initial)

    #7. Make sure that all land cells are filled
    subcat5 = pcr.cover(subcat4,land)
    
    #8. Each water region is given a unique ID number
    waterregion0 = pcr.nominal(subcat5)
    waterregion = pcr.clump(waterregion0)


    #9. Save the water region map
    pcr.report(waterregion,output_wr)
    
    return waterregion, subcat1
                          
def main_script():
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    main_script()                         

