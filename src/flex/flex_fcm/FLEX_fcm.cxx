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

#include "FLEX_fcm.h"

#include <vector>
#include <map>
#include <iostream>
#include <numeric>
#include <math.h>
#include <time.h>
#include <cmath>
#include <stdlib.h>


FCM::FCM( std::vector< double >& data, int cluster_n,
          double expo, int max_iter, double min_improv ) 
	: data(data), 
		cluster_n(cluster_n),
		expo(expo),
		max_iter(max_iter),
		min_improv(min_improv),
		data_n( data.size() )
{		
	
	U = CreateVector();
	Unew = CreateVector();

	center.resize( cluster_n );

	obj_fcn = 0;
			
	InitU();
}

FCM::~FCM()
{
	DeleteVector( U );
	DeleteVector( Unew );
}

void FCM::InitU()
{
	// Create partition matrix.  
		// Number of columns equal to number of datapoints to be clustered.
		// Number of rows equal to the number of clusters.
		// Thus, each column holds all the fuzzy membership values
		// for an individual voxel.

	// Randomly initialize partition matrix ( 1 <= U > 0 ).
  srand( time(NULL) ); 
  for ( int i = 0; i < cluster_n; ++i )
		for ( int j = 0; j < data_n; ++j )
			U[i]->operator [](j) = ( rand() % 100001 / 100000.0 ) + 0.0005;
			//** add 0.0005 to ensure no col in U is all 0s **//
	
	// Normalize partition matrix such that the  
	// membership values in each column sum to unity.
	for ( int j = 0; j < data_n; ++j )
	{
		double sum = 0;
		for ( int i =0; i < cluster_n; ++i )
			sum += U[i]->operator [](j);

		for ( int i = 0; i < cluster_n; ++i )
			 U[i]->operator [](j) /= sum ; 


  } // for 
    

} // InitU(...)


void FCM::StepFcm()
{
	std::vector< std::vector< double >* > mf = CreateVector();
	
	// mf = U.^expo
	// Exponential modification of partition matrix.
	for ( int i = 0; i < cluster_n; ++i )
  {
    for ( int j = 0; j < data_n; ++j )
		{	
			mf[i]->operator [](j) = pow( U[i]->operator [](j), 2 );
		} 
  } 
		
	// center = mf*data./((ones(size(data, 2), 1)*sum(mf'))');  
  for ( int i = 0; i < cluster_n; ++i ) 
	{
    center[i] = inner_product( mf[i]->begin(), mf[i]->end(), data.begin(), 0.0 ) / accumulate( mf[i]->begin(), mf[i]->end(), 0.0 );
	}
  	
  
	// dist = distfcm(center, data); 
	std::vector< std::vector< double >* > dist = CreateVector();
	for (int i = 0; i < cluster_n; ++i )
  {
    for (int j = 0; j < data_n; ++j )
		{
			dist[i]->operator [](j) = fabs( data[j] - center[i] );
		}
  }


	// obj_fcn = sum(sum((dist.^2).*mf));		
	//  Update objective function.
  std::vector< std::vector< double >* > temp = CreateVector();
	for ( int i = 0; i < cluster_n; ++i )
  {
		for ( int j = 0; j < data_n; ++j )
		{
			temp[i]->operator [](j) = ( dist[i]->operator [](j) * dist[i]->operator [](j) ) / (mf[i]->operator [](j)+.1);
		}
  }

	obj_fcn = 0.0;
	for ( int i = 0; i < cluster_n; ++i )
	  obj_fcn = obj_fcn + accumulate( temp[i]->begin(), temp[i]->end(), 0.0 );

  
	// tmp = dist.^(-2/(expo-1));
	for ( int i = 0; i < cluster_n; ++i )
  {	
    for ( int j = 0; j < data_n; ++j )
		{
			temp[i]->operator [](j) = pow ( dist[i]->operator [](j), ( -2 / (expo-1) ) );
		}
  }

	// U_new = tmp./(ones(cluster_n, 1)*sum(tmp));
	for ( int j = 0; j < data_n; ++j )
	{
		double sum = 0;
		for ( int i = 0; i < cluster_n; ++i )
			sum = sum + temp[i]->operator [](j);
	
		for ( int i = 0; i < cluster_n; ++i )
		{
			Unew[i]->operator [](j) = temp[i]->operator [](j) / sum;
		}
	}

	DeleteVector( mf );
	DeleteVector( dist );
	DeleteVector( temp );
}

void FCM::RunFcm()
{
	double prev_obj_fcn;

	for ( int i = 0; i < max_iter; ++i )
	{
    
    std::cout<<std::left<<std::fixed;
    
    std::cout << "Iteration count = " << i + 1 << ", Objective function = " << obj_fcn << std::endl;
		
    prev_obj_fcn = obj_fcn;

		StepFcm();
		
    if ( i != 0 )
		  if ( fabs( obj_fcn - prev_obj_fcn ) < min_improv )
			  break; 

		for ( int i = 0; i < cluster_n; ++i )
			copy( Unew[i]->begin(), Unew[i]->end(), U[i]->begin() );

	}
}

std::vector< std::vector< double >* >  FCM::DistFcm()
{
	std::vector< std::vector< double >* > out = CreateVector();
	
	for (int i = 0; i < cluster_n; ++i )
		for (int j = 0; j < data_n; ++j )
			out[i]->operator [](j) = std::abs( data[j] - center[i] );

	return out;		
}
	
std::vector< std::vector< double >* > FCM::CreateVector()
{
	std::vector< std::vector< double >* > newVector( cluster_n );

	for ( int i = 0; i < cluster_n; ++i )
	{
		newVector[i] = new std::vector< double >( data_n );
	}
	return newVector;
}

void FCM::DeleteVector( std::vector< std::vector< double >* > oldVector )
{
	for ( int i = 0; i < cluster_n; ++i )
		delete oldVector[i];			
}

 std::vector< double>*  FCM::GetU(int classNum)
{
	typedef std::map< double, int > DoubleIntMap;

	DoubleIntMap map;

	for ( int i = 0; i < cluster_n; ++ i )
		map[ center[i] ] = i;

	DoubleIntMap::iterator it;
	it = map.begin();
	
	for ( int i = 0; i < classNum - 1; ++it, ++i  )
		;

	
	return Unew[ it->second ];
}
