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

#include "egRunningStats.h"

namespace eg
{
  
  void RunningStats::Push( double x )
  {
    ++m_n;

    if ( m_n == 1 )
    {
      m_oldM = m_newM = x;
      m_oldV = 0.0;
    }
    else
    {
      m_newM = m_oldM + ( x - m_oldM ) / m_n;
      m_newV = m_oldV + ( x - m_oldM ) * ( x - m_newM );

      // setup for next iteration
      m_oldM = m_newM;
      m_oldV = m_newV;
    }

  };

}