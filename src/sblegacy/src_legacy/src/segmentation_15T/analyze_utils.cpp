/*=========================================================================
  
  Copyright 2002, 2006-2010, Dr. Sandra Black
  Linda C. Campbell Cognitive Neurology Unit
  Sunnybrook Health Sciences Center
  
  This file is part of the Sunnybrook Image Software Processing (SIPS) package

  SIPS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#include <iostream.h>
#include <fstream.h>
#include <vector.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "dbh.h"


void make_hdr(char *hdr_name,
              short x_size, short y_size, short z_size,
              int num, char *datatype,
              int max, int min,
              float x_pixsize, float y_pixsize, float z_pixsize)
    {
    int i;
    struct dsr hdr;
    FILE *fp;
    static char DataTypes[9][12] = {"UNKNOWN", "BINARY", "CHAR", "SHORT", "INT",
                                    "FLOAT", "COMPLEX", "DOUBLE", "RGB"};
                                                            
    static int DataTypeSizes[9] = {0,1,8,16,32,32,64,64,24};
    
    memset(&hdr,0, sizeof(struct dsr));
    for(i=0;i<8;i++)
        hdr.dime.pixdim[i]=0.0;
   
    hdr.dime.vox_offset = 0.0;
    hdr.dime.roi_scale   = 1.0;
    hdr.dime.funused1    = 0.0;
    hdr.dime.funused2    = 0.0;
    hdr.dime.cal_max     = 0.0;
    hdr.dime.cal_min     = 0.0;
  
   
    hdr.dime.datatype = -1;
    char temp_datatype[10];
    strcpy(temp_datatype, datatype);
    for(i=1;i<=8;i++)
        if(!strcmp(temp_datatype,DataTypes[i]))
        {
                hdr.dime.datatype = (1<<(i-1));
                hdr.dime.bitpix = DataTypeSizes[i];
                break;
        }
 
    if((fp=fopen(hdr_name,"w"))==0)
        {
        printf("unable to create: %s\n", hdr_name);
        exit(0);
    }

    hdr.dime.dim[0] = 4;  /* all Analyze images are taken as 4 dimensional */
    hdr.hk.regular = 'r';
    hdr.hk.sizeof_hdr = sizeof(struct dsr);

    hdr.dime.dim[1] = x_size;  /* slice width  in pixels */
    hdr.dime.dim[2] = y_size;  /* slice height in pixels */
    hdr.dime.dim[3] = z_size;  /* volume depth in slices */
    hdr.dime.dim[4] = num;  /* number of volumes per file */

    hdr.dime.glmax  = max;  /* maximum voxel value  */
    hdr.dime.glmin  = min;  /* minimum voxel value */
{
    hdr.dime.pixdim[1] = x_pixsize; /* voxel width in mm per pixel */

        hdr.dime.pixdim[2] = y_pixsize; /* voxel height in mm per pixel */
        hdr.dime.pixdim[3] = z_pixsize; /* voxel depth in mm per pixel
        */
}
/*   Assume zero offset in .img file, byte at which pixel
       data starts in the image file */

    hdr.dime.vox_offset = 0.0; 
    
/*   Planar Orientation;    */
/*   Movie flag OFF: 0 = transverse, 1 = coronal, 2 = sagittal
     Movie flag ON:  3 = transverse, 4 = coronal, 5 = sagittal  */  

    hdr.hist.orient     = 0;  
    
/*   up to 3 characters for the voxels units label; i.e. 
        mm., um., cm.               */

    strcpy(hdr.dime.vox_units," ");
   
/*   up to 7 characters for the calibration units label; i.e. HU */

    strcpy(hdr.dime.cal_units," ");  
    
/*     Calibration maximum and minimum values;  
       values of 0.0 for both fields imply that no 
    calibration max and min values are used    */

    hdr.dime.cal_max = 0.0; 
    hdr.dime.cal_min = 0.0;

    fwrite(&hdr,sizeof(struct dsr),1,fp);
    fclose(fp);
}
////////////////////////////////////////////////////////////////////////////
