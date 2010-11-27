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

#ifndef egRunningStats_h__
#define egRunningStats_h__

#include <cmath>

namespace eg
{

  class RunningStats
  {
  public:
    RunningStats() : m_n( 0 ) {};

    void Clear() { m_n = 0; }
    void Push( double x );

    int NumDataValues() const        { return m_n; }
    double Mean() const              { return ( m_n > 0 ) ? m_newM : 0.0; }
    double Variance() const          { return ( m_n > 1 ) ? m_newV / ( m_n-1 ) : 0.0; }
    double StandardDeviation() const { return std::sqrt( this->Variance() ); }

  private:
    int    m_n;
    double m_oldM, m_newM;
    double m_oldV, m_newV;
  };

}
#endif