from cc3d.core.PySteppables import *
from pathlib import Path
from random import uniform, seed, randint, sample
from random import randrange
from numpy.random import normal
from numpy import mean


import csv
import numpy as np
import math
import random



class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):

        """
        any code in the start function runs before MCS=0
        """
        # inititiaize a single cell at the initial simulation of MCS=0
        
        self.cell_field[150:158, 150:158, 0] = self.new_cell(self.NETWORKS)
        #self.cell_field[150:158, 150:158, 0] = self.new_cell(self.SPHEROIDS)
   
        
        # set cell target volume and lambda volume
        for cell in self.cell_list:
            cell.targetVolume = 64
            cell.lambdaVolume = 2.0
        

        ## wrtie data to .csv files        
        output_dir = self.output_dir
        header1 = ['mcs','cellType','x','y']
        header2 = ['mcs','cellid','cellType','x','y']
        header3 = ['mcs','circularity','cell_invasion']
        
        if output_dir is not None:
            output_path1 = Path(output_dir).joinpath('cellposition.csv')
            output_path2 = Path(output_dir).joinpath('cell_typeposition.csv')
            output_path3 = Path(output_dir).joinpath('celldata.csv')
            
            with open(output_path1,'a',encoding='UTF8',newline='') as f1:
                writer1 = csv.writer(f1)
                writer1.writerow(header1)
                
            with open(output_path2,'a',encoding='UTF8',newline='') as f2:
                writer2 = csv.writer(f2)
                writer2.writerow(header2)
        
            with open(output_path3,'a',encoding='UTF8',newline='') as f3:
                writer1 = csv.writer(f3)
                writer1.writerow(header3)
        
        # initialize the plot for the cell circularity
        self.plot_win1 = self.add_new_plot_window(title='Circularity',
                                                 x_axis_title='Time (days)',
                                                 y_axis_title='Circularity', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)

        self.plot_win1.add_plot("Circularity", style='Dots', color='red', size=1)
       
        # initialize the plot for the cell Invasion 
        self.plot_win = self.add_new_plot_window(title='Invasion',
                                                 x_axis_title='Time(days)',
                                                 y_axis_title='Distance (micron)', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)

        self.plot_win.add_plot("Invasion", style='Dots', color='red', size=1)
        

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        cell_count = len(self.cell_list)
        for cell in self.cellList:
            if cell.volume<=128 and cell_count==1:
                cell.targetVolume = cell.targetVolume + 64./800
                cell.lambdaVolume = 2.0
            if cell.volume<=128 and cell_count>1:
                cell.targetVolume = cell.targetVolume + 64./2000
                cell.lambdaVolume = 2.0
            #print(cell.targetVolume, cell.type)
        

        distance_max = 0
        corrected_length = 0
        for cell1 in self.cellList:
            for cell2 in self.cellList:               
                distance = np.sqrt((cell1.xCOM-cell2.xCOM)**2 + (cell1.yCOM-cell2.yCOM)**2)
                if distance > distance_max:
                    distance_max = distance
                    # adding the radius of two cell with center to center distance
                    corrected_length = np.sqrt(cell1.volume/3.1416) + np.sqrt(cell2.volume/3.1416)
                if corrected_length == 0:
                    corrected_length = np.sqrt(cell1.volume/3.1416) + np.sqrt(cell2.volume/3.1416)
        
        total_volume = 0
        for cell in self.cellList:
            total_volume = total_volume + cell.volume
        
        
        #tot_v1 = (distance_max*distance_max*3.1416)/4
        
        corrected_circular_volume = ((distance_max+corrected_length)*(distance_max+corrected_length)*3.1416)/4
        
        distance_Invasion = 0
        for cell_migrated in self.cellList:
            distance1 = np.sqrt((155-cell_migrated.xCOM)*(155-cell_migrated.xCOM) + (155-cell_migrated.yCOM)*(155-cell_migrated.yCOM))
            if distance1 > distance_Invasion:
                distance_Invasion = distance1
                
        circularity = total_volume/corrected_circular_volume
        cell_invasion = distance_Invasion+(corrected_length)/2
        
        # print(distance_max, corrected_length, total_volume)
        # #print(total_volume/corrected_circular_volume, distance_Invasion)
        # print(circularity, distance_Invasion)
        
        
       
        self.plot_win.add_data_point("Invasion", mcs/1000, distance_Invasion*2)
        self.plot_win1.add_data_point("Circularity", mcs/1000, circularity)
      
    
                
        ## write data to .csv files 
        output_dir = self.output_dir
    
        
        if output_dir is not None:
            output_path1 = Path(output_dir).joinpath('cellposition.csv')
            output_path2 = Path(output_dir).joinpath('cell_typeposition.csv')
            output_path3 = Path(output_dir).joinpath('celldata.csv')
             
            with open(output_path1,'a',encoding='UTF8',newline='') as f1:
                writer1 = csv.writer(f1)
                writer1.writerow([mcs,cell.type,cell.xCOM,cell.yCOM])
                
            with open(output_path2,'a',encoding='UTF8',newline='') as f2:
                writer2 = csv.writer(f2)
                for cell in self.cell_list:
                    writer2.writerow([mcs,cell.id,cell.type,cell.xCOM,cell.yCOM])

            with open(output_path3,'a',encoding='UTF8',newline='') as f3:
                writer3 = csv.writer(f3)
                writer3.writerow([mcs,(total_volume/corrected_circular_volume),(distance_Invasion)])   #+(corrected_length)/2
                      

               
    def finish(self):
        # save plots
    
                
        if self.output_dir is not None:
            png_output_path11 = Path(self.output_dir).joinpath("Invasion.png")
            # here we SPHEROIDSNETWORKSify size of the image saved - default is 400 x 400
            self.plot_win.save_plot_as_png(png_output_path11, 1000, 1000) 
            
            
            png_output_path12 = Path(self.output_dir).joinpath("Circularity.png")
            # here we SPHEROIDSNETWORKSify size of the image saved - default is 400 x 400
            self.plot_win1.save_plot_as_png(png_output_path12, 1000, 1000)
        
class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
    
        for cell in self.cell_list:
            cell.targetVolume += 1        

        # # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        

        # field = self.field.CHEMICAL_FIELD_NAME
        
        # for cell in self.cell_list:
            # concentrationAtCOM = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]

            # # you can use here any fcn of concentrationAtCOM
            # cell.targetVolume += 0.01 * concentrationAtCOM       

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
        self.set_parent_child_position_flag(-1)

    def step(self, mcs):
        cells_to_divide=[]  
        for cell_1 in self.cell_list:
            if cell_1.volume>=128:
                cells_to_divide.append(cell_1)
        for cell_1 in cells_to_divide:
            self.divide_cell_random_orientation(cell_1)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        #self.parent_cell.targetVolume /= 2.0
        #self.lengthConstraintLocalFlexPlugin.setLengthConstraintData(self.parent_cell, 1.4, 30)

        self.clone_parent_2_child()
        self.childCell.targetVolume = 64
        self.childCell.lambdaVolume = 2
        self.parentCell.targetVolume = 64
        self.parentCell.lambdaVolume = 2

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        value1 = random.uniform(0, 1)
        value2 = random.uniform(0, 1)
        
        if self.parent_cell.type==1 and value1>=0.5:
            self.child_cell.type=1
        if self.parent_cell.type==1 and value1<0.5:
            self.child_cell.type=1
        
        if self.parent_cell.type==2 and value2>=0.5:
            self.child_cell.type=2
        if self.parent_cell.type==2 and value2<0.5:
            self.child_cell.type=2

        