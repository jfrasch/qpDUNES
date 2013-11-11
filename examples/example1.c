/*
 *	This file is part of qpDUNES.
 *
 *	qpDUNES -- A DUal NEwton Strategy for convex quadratic programming.
 *	Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau et al. 
 *	All rights reserved.
 *
 *	qpDUNES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpDUNES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpDUNES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file examples/example1.c
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 *
 *	Very simple example for testing qpDUNES.
 */



#include <qpDUNES.h>

#define INFTY 1.0e12

int main( )
{
	int i;
	boolean_t isLTI;
	
	
	/* set dimensions */
	unsigned int nI = 2;		/* number of stages */
	unsigned int nX = 3;		/* number of states */
	unsigned int nU = 2;		/* number of controls */
	unsigned int* nD = 0;		/* number of affine constraints */
	
	
	/* specify problem data */
	double Q[3*3] =
		{	1.0, 0.0, 0.0, 
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0	};			
	double R[2*2] =
		{	1.0, 0.0,
			0.0, 1.0	};
	double *S=0;
	
	double* P = Q;
	
	double A[3*3] =
		{	1.0, 0.0, 0.0, 
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0	};
	double B[3*2] =
		{	1.0, 0.0,
			0.0, 1.0,
			1.0, 1.0	};
	double c[3] = 
		{	5.0,
			5.0,
			5.0	};
	
	double xLow[3] = { -INFTY, -INFTY, -INFTY };
	double xUpp[3] = {  INFTY,  INFTY,  INFTY };
	double uLow[2] = { -INFTY, -INFTY };
	double uUpp[2] = {  INFTY,  INFTY };
	
	
	/* qpData struct */
	qpData_t qpData;


	/* memory allocation */
	qpDUNES_setup( &qpData, nI, nX, nU, nD, 0 );	/* passing 0 in the last argument sets the default QP options */


	/* manual setup of intervals */
	for( i=0; i<nI; ++i )
	{
		qpDUNES_setupSimpleBoundedInterval(  &qpData, qpData.intervals[i],Q,R,S, A,B,c, xLow,xUpp,uLow,uUpp );
	}
	qpDUNES_setupSimpleBoundedInterval(  &qpData, qpData.intervals[nI], P,0,0, 0,0,0, xLow,xUpp,0,0 );

	qpDUNES_setupAllLocalQPs( &qpData, isLTI=QPDUNES_TRUE );	/* determine local QP solvers and set up auxiliary data */


	/* solve problem */
	qpDUNES_solve( &qpData );
	
	
	/* write out solution */
	for( i=0; i<nI; ++i )
	{
		qpDUNES_printMatrixData( qpData.intervals[i]->z.data, 1, nX+nU, "z[%d]:", i );
	}
	qpDUNES_printMatrixData( qpData.intervals[nI]->z.data, 1, nX, "z[%d]:", i );
	
	qpDUNES_cleanup( &qpData );
	
	printf( "example1 done.\n" );
	
	return 0;
}


/*
 *	end of file
 */
