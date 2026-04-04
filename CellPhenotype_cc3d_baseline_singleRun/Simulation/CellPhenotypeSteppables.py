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

cell_type = "Spheroid" # "Spheroid", "Network"
 
class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        self.original_cell = None
        self.original_cell_volume = None
      
      

    def start(self):

        """
        any code in the start function runs before MCS=0
        """
        # inititiaize a single cell at the initial simulation of MCS=0
        
        selected_cell_type = "SPHEROID"  # or "NETWORK"
        #selected_cell_type = "NETWORK"

        # Cell width logic
        cell_type_map = {
            "SPHEROID": (self.SPHEROIDS, 9),
            "NETWORK": (self.NETWORKS, 9)
        }

        cellType, cellWidth = cell_type_map[selected_cell_type]

        xMid = self.dim.x / 2
        yMid = self.dim.y / 2

        newCell = self.new_cell(cellType)
        self.original_cell = newCell  # Store reference to the original cell
    

        xCellCenter = int(xMid)
        yCellCenter = int(yMid)

        for ix in range(xCellCenter - cellWidth, xCellCenter + cellWidth + 1):
            for iy in range(yCellCenter - cellWidth, yCellCenter + cellWidth + 1):
                if math.sqrt((ix - xCellCenter) ** 2 + (iy - yCellCenter) ** 2) <= cellWidth / 2:
                    #self.cell_field[ix:ix + 1, iy:iy + 1, 0] = newCell
                    self.cell_field[int(ix), int(iy), 0] = newCell

        newCell.targetVolume = newCell.volume
        newCell.lambdaVolume = 2.0

        #print(f"Initialized cell of type {selected_cell_type}, target volume = {newCell.targetVolume}")
        #print(f"Initialized cell of type {selected_cell_type}, X-axis center = {self.original_cell.xCOM}")
        #print(f"Initialized cell of type {selected_cell_type}, Y-axis center = {self.original_cell.yCOM}")



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
        firstDoubleTime = 700
        secondDoubleTime = 3000 # 3000 for spheroid,for Network

        # Update volume for growth
        cell_count = len(self.cell_list)
        for cell in self.cellList:
            if cell.volume <= self.original_cell.volume * 2 and cell_count == 1:
                cell.targetVolume += self.original_cell.volume / firstDoubleTime
                cell.lambdaVolume = 2.0
            elif cell.volume <= self.original_cell.volume * 2 and cell_count > 1:
                cell.targetVolume += self.original_cell.volume / secondDoubleTime
                cell.lambdaVolume = 2.0

        # Circularity calculations
        distance_max = 0
        corrected_length = 0
        for cell1 in self.cellList:
            for cell2 in self.cellList:
                distance = np.sqrt((cell1.xCOM - cell2.xCOM) ** 2 + (cell1.yCOM - cell2.yCOM) ** 2) 
                if distance > distance_max:
                    distance_max = distance
                    corrected_length = np.sqrt(cell1.volume / np.pi) + np.sqrt(cell2.volume / np.pi)
                if corrected_length == 0:
                    corrected_length = np.sqrt(cell1.volume / np.pi) + np.sqrt(cell2.volume / np.pi)

        total_volume = sum(cell.volume for cell in self.cellList)
        corrected_circular_volume = ((distance_max + corrected_length) ** 2 * np.pi) / 4

        # Invasion calculation
        distance_Invasion = 0
        if self.original_cell is not None:
            for cell_migrated in self.cellList:
                dx = cell_migrated.xCOM - self.original_cell.xCOM
                dy = cell_migrated.xCOM - self.original_cell.yCOM #force initial cell_migrated.yCOM = cell_migrated.xCOM
                distance = np.sqrt(dx ** 2 + dy ** 2)  
                if distance > distance_Invasion:
                    distance_Invasion = distance #- np.sqrt(cell_migrated.volume/np.pi)

        circularity = total_volume / corrected_circular_volume
        cell_invasion = distance_Invasion # + (corrected_length / 2)
        
       
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
           

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
        self.set_parent_child_position_flag(-1)
        self.initialCellVolume = 69.0 #63.64
        
  
    def step(self, mcs):
        initialCellVolume = 69.0 #63.64
        cells_to_divide=[]  
        for cell_1 in self.cell_list:
            if cell_1.volume>= initialCellVolume *2:  #initial cell volume times two
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
        self.childCell.targetVolume = self.initialCellVolume #64,51 initial cell volume
        self.childCell.lambdaVolume = 2
        self.parentCell.targetVolume = self.initialCellVolume
        self.parentCell.lambdaVolume = 2

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
                
  
        
        if self.parent_cell.type==2 and cell_type=="Network":
            self.child_cell.type=2
            self.parent_cell.type=2
        
        if self.parent_cell.type==1 and cell_type=="Spheroid":
            self.child_cell.type=1
            self.parent_cell.type=1
        
        
        
        
        # value1 = random.uniform(0, 1)
        # value2 = random.uniform(0, 1)
        
        # if self.parent_cell.type==1 and value1>=0.5:
            # self.child_cell.type=1
        # if self.parent_cell.type==1 and value1<0.5:
            # self.child_cell.type=1
        
        # if self.parent_cell.type==2 and value2>=0.5:
            # self.child_cell.type=2
        # if self.parent_cell.type==2 and value2<0.5:
            # self.child_cell.type=2
