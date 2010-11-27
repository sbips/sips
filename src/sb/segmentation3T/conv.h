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

#ifndef CONV_H
#define CONV_H
/* converts integer */

unsigned short suiconv( unsigned short iOld )
{
        char* pcOrig;
        char* pcNew;
        unsigned short iNew;

        pcOrig = (char *) &iOld;
        pcNew = (char *) &iNew;

        pcNew[0] = pcOrig[1];
        pcNew[1] = pcOrig[0];

        return iNew;
}

int iconv( int iOld )
{
        char* pcOrig;
        char* pcNew;
        int i, j, iNew;

        pcOrig = (char *) &iOld;
        pcNew = (char *) &iNew;

        j = sizeof(int)-1;
        for (i = 0; i<sizeof(int); i++) {
                pcNew[j] = pcOrig[i];
                j--;
                }

        return iNew;
}

int lconv( long iOld )
{
        char* pcOrig;
        char* pcNew;
        int i, j;
	long iNew;

        pcOrig = (char *) &iOld;
        pcNew = (char *) &iNew;

        j = sizeof(long)-1;
        for (i = 0; i<sizeof(long); i++) {
                pcNew[j] = pcOrig[i];
                j--;
                }

        return iNew;
}

/* converts short integer */
short int siconv( short int iOld )
{
        char* pcOrig;
        char* pcNew;
        short int iNew;

        pcOrig = (char *) &iOld;
        pcNew = (char *) &iNew;

        pcNew[0] = pcOrig[1];
        pcNew[1] = pcOrig[0];

        return iNew;
}

/* converts float */
float fconv( float fOld )
{
        char* pcOrig;
        char* pcNew;
        int i, j;
        float fNew;

        pcOrig = (char *) &fOld;
        pcNew = (char *) &fNew;

        j = sizeof(float)-1;
        for (i = 0; i<sizeof(float); i++) {
                pcNew[j] = pcOrig[i];
                j--;
                }

        return fNew;
}

#endif

