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

//Last Modified: 10/12/2008

#ifndef _FCM_H_
#define _FCM_H_

#include <vector>


class FCM
{

public:
	
	FCM( std::vector< double >& data, int cluster_n, double expo, int max_iter, double min_improv ); 

	~FCM();
	
	void RunFcm();

	std::vector< double >* GetU(int classNum);

	
private:

	void InitU();

	void StepFcm();

	std::vector< std::vector< double >* >  DistFcm();
	
	std::vector< std::vector< double >* > CreateVector();

	void DeleteVector( std::vector< std::vector< double >* > oldVector );


	std::vector< double > const data;
	
  const int    cluster_n;
	const double expo;
	const int    max_iter;
	const double min_improv;
	const int    data_n;

	double obj_fcn;

	std::vector< std::vector< double >* > U;
	std::vector< std::vector< double >* > Unew;

	std::vector< double > center;

};


#endif
