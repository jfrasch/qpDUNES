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
 *	\file interfaces/matlab/qpDUNES.cpp
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 *
 *	Interface for Matlab(R) that enables to call qpDUNES as a MEX function.
 *
 */


extern "C" {
#include <qpDUNES.h>
}
#include <qpDUNES_matlab_utils.cpp>

#include "mex.h"
#include "matrix.h"
#include "string.h"


/* global pointer to qpDUNES object */
static qpData_t* qpDataGlobal = 0;



/*
 *	q p D U N E S _ s a f e C l e a n u p M a t l a b
 */
static void qpDUNES_safeCleanupMatlab(
										 )
{
	if (qpDataGlobal != 0) {
		qpDUNES_cleanup( qpDataGlobal );
		delete qpDataGlobal;
		qpDataGlobal = 0;
	}
}
/*<<< END OF qpDUNES_cleanupMatlab */



/*
 *	q p D U N E S _ s e t u p M a t l a b
 */
void qpDUNES_setupMatlab( qpData_t** qpDataPtr,
						  uint_t nI,
						  uint_t nX,
						  uint_t nU,
						  uint_t* nD,
						  const mxArray* const optionsPtr
						  )
{
	qpDUNES_safeCleanupMatlab( );
	*qpDataPtr = new qpData_t;
	qpOptions_t qpOptions = qpDUNES_setupDefaultOptions();

	/* set up user options */
    if ( optionsPtr != 0 ) {
    	qpDUNES_setupOptionsMatlab( &(qpOptions), optionsPtr );
    }

	/* allocate qpDUNES data */
	qpDUNES_setup( *qpDataPtr, nI, nX, nU, nD, &(qpOptions) );

	/* register qpDUNES memory for MATLAB-triggered clearing */
	mexAtExit(qpDUNES_safeCleanupMatlab);
}
/*<<< END OF qpDUNES_setupMatlab */



/*
 *	i n i t Q P
 */
void initQP( int nrhs, const mxArray* const prhs[] )
{
	/* 0) VARIABLE DECLARATIONS: */
	return_t statusFlag;

	mxArray* nIuserPtr = 0;
	real_t* nIuser = 0;

	mxArray* Hptr = 0;
	real_t* H_in = 0;
	mxArray* Pptr = 0;
	real_t* P = 0;
	mxArray* gPtr = 0;
	real_t* g = 0;

	mxArray* Cptr = 0;
	real_t* C_F = 0;
	mxArray* cPtr = 0;
	real_t* c = 0;

	mxArray* zLowPtr = 0;
	real_t* zLow = 0;
	mxArray* zUppPtr = 0;
	real_t* zUpp = 0;

	mxArray* Dptr = 0;
	real_t* D_F = 0;
	mxArray* dLowPtr = 0;
	real_t* dLow = 0;
	mxArray* dUppPtr = 0;
	real_t* dUpp = 0;

	const mxArray* optionsPtr = 0;


	/* I) GET OUTPUT POINTERS: */
	nIuserPtr = (mxArray*)prhs[0];
	Hptr = (mxArray*)prhs[1];	/* to get/check dimensions */
	Pptr = (mxArray*)prhs[2];
	gPtr = (mxArray*)prhs[3];
	Cptr = (mxArray*)prhs[4];
	cPtr = (mxArray*)prhs[5];
	zLowPtr = (mxArray*)prhs[6];
	zUppPtr = (mxArray*)prhs[7];
	Dptr = (mxArray*)prhs[8];
	dLowPtr = (mxArray*)prhs[9];
	dUppPtr = (mxArray*)prhs[10];
	if ( ( nrhs == 12 ) &&					/* check whether options are specified */
		 ( !mxIsEmpty(prhs[nrhs-1]) ) &&
		 ( mxIsStruct(prhs[nrhs-1]) ) )
	{
		optionsPtr = prhs[nrhs-1];
	}
	/* exatract data */
	nIuser = (real_t*)  mxGetPr( nIuserPtr );
	H_in = (real_t*) mxGetPr( Hptr );
	P = (real_t*) mxGetPr( Pptr );
	g = (real_t*) mxGetPr( gPtr );
	C_F = (real_t*) mxGetPr( Cptr );
	c = (real_t*) mxGetPr( cPtr );
	zLow = (real_t*) mxGetPr( zLowPtr );
	zUpp = (real_t*) mxGetPr( zUppPtr );
	D_F = (real_t*) mxGetPr( Dptr );
	dLow = (real_t*) mxGetPr( dLowPtr );
	dUpp = (real_t*) mxGetPr( dUppPtr );


	/* 2) check pointer consistency */
	if ( nIuser == 0 )	mexErrMsgTxt( "ERROR (qpDUNES): Number of Intervals nI undefined!" );
	if ( mxGetM( nIuserPtr )*mxGetN( nIuserPtr ) != 1 )	mexErrMsgTxt( "ERROR (qpDUNES): Number of Intervals nI wrongly specified!" );
	if ( Cptr == 0 )	mexErrMsgTxt( "ERROR (qpDUNES): C matrix undefined!" );


	/* 3) get dimensions */
	uint_t nI = (uint_t)(*nIuser);
	uint_t nColsH, nColsC;
	uint_t nX, nU, nZ;
	uint_t* nD = new uint_t[nI+1];

	nZ = (uint_t) mxGetM( Hptr );			/* get number of rows */
	nColsH = (uint_t) mxGetN( Hptr );		/* get number of columns */

	nX = (uint_t) mxGetM( Cptr );		/* get number of rows */
	nColsC = (uint_t) mxGetN( Cptr );	/* get number of columns */

	nU = nZ - nX;

	/* TODO: deal with affine constraints! */
	for ( uint_t ii=0; ii<nI; ++ii ) {
		/* 		nD[ii] = nZ; */
		nD[ii] = 0;
	}
	/* 	nD[nI] = nX;	*/
	nD[nI] = 0;


	/* 4) check dimension consistency */
	if ( nColsH != nI*nZ )	{
		mexPrintf( "[qpDUNES] Dimensions of data received:\n  nI  : %d\n  H   : %d x %d\n  P   : %d x %d\n  g   : %d x %d\n  C   : %d x %d\n  c   : %d x %d\n  zLow: %d x %d\n  zUpp: %d x %d\n",
									nI,
									mxGetM( Hptr ), mxGetN( Hptr ),
									mxGetM( Pptr ), mxGetN( Pptr ),
									mxGetM( gPtr ), mxGetN( gPtr ),
									mxGetM( Cptr ), mxGetN( Cptr ),
									mxGetM( cPtr ), mxGetN( cPtr ),
									mxGetM( zLowPtr ), mxGetN( zLowPtr ),
									mxGetM( zUppPtr ), mxGetN( zUppPtr )
									);
		mexErrMsgTxt( "[qpDUNES] Error: Dimensions of H are inconsistent." );
	}
	if ( ( P != 0 ) && ( ( mxGetM( Pptr ) != nX ) || ( mxGetN( Pptr ) != nX ) ) ) {
		mexPrintf( "[qpDUNES] Dimensions of data received:\n  nI  : %d\n  H   : %d x %d\n  P   : %d x %d\n  g   : %d x %d\n  C   : %d x %d\n  c   : %d x %d\n  zLow: %d x %d\n  zUpp: %d x %d\n",
									nI,
									mxGetM( Hptr ), mxGetN( Hptr ),
									mxGetM( Pptr ), mxGetN( Pptr ),
									mxGetM( gPtr ), mxGetN( gPtr ),
									mxGetM( Cptr ), mxGetN( Cptr ),
									mxGetM( cPtr ), mxGetN( cPtr ),
									mxGetM( zLowPtr ), mxGetN( zLowPtr ),
									mxGetM( zUppPtr ), mxGetN( zUppPtr )
									);
		mexErrMsgTxt( "[qpDUNES] Error: Detected P, but dimensions are inconsistent." );
	}
	if ( ( g != 0 ) && ( mxGetM( gPtr )*mxGetN( gPtr ) != nI*nZ + nX ) ) {
		mexPrintf( "[qpDUNES] Dimensions of data received:\n  nI  : %d\n  H   : %d x %d\n  P   : %d x %d\n  g   : %d x %d\n  C   : %d x %d\n  c   : %d x %d\n  zLow: %d x %d\n  zUpp: %d x %d\n",
									nI,
									mxGetM( Hptr ), mxGetN( Hptr ),
									mxGetM( Pptr ), mxGetN( Pptr ),
									mxGetM( gPtr ), mxGetN( gPtr ),
									mxGetM( Cptr ), mxGetN( Cptr ),
									mxGetM( cPtr ), mxGetN( cPtr ),
									mxGetM( zLowPtr ), mxGetN( zLowPtr ),
									mxGetM( zUppPtr ), mxGetN( zUppPtr )
									);
		mexErrMsgTxt( "[qpDUNES] Error: Detected g, but dimensions are inconsistent." );
	}
	if ( nColsC != nI*nZ ) {
		mexPrintf( "[qpDUNES] Dimensions of data received:\n  nI  : %d\n  H   : %d x %d\n  P   : %d x %d\n  g   : %d x %d\n  C   : %d x %d\n  c   : %d x %d\n  zLow: %d x %d\n  zUpp: %d x %d\n",
									nI,
									mxGetM( Hptr ), mxGetN( Hptr ),
									mxGetM( Pptr ), mxGetN( Pptr ),
									mxGetM( gPtr ), mxGetN( gPtr ),
									mxGetM( Cptr ), mxGetN( Cptr ),
									mxGetM( cPtr ), mxGetN( cPtr ),
									mxGetM( zLowPtr ), mxGetN( zLowPtr ),
									mxGetM( zUppPtr ), mxGetN( zUppPtr )
									);
		mexErrMsgTxt( "[qpDUNES] Error: Dimensions of C are inconsistent." );
	}
	if ( ( c != 0 ) && ( mxGetM( cPtr )*mxGetN( cPtr ) != nI*nX ) ) {
		mexPrintf( "[qpDUNES] Dimensions of data received:\n  nI  : %d\n  H   : %d x %d\n  P   : %d x %d\n  g   : %d x %d\n  C   : %d x %d\n  c   : %d x %d\n  zLow: %d x %d\n  zUpp: %d x %d\n",
									nI,
									mxGetM( Hptr ), mxGetN( Hptr ),
									mxGetM( Pptr ), mxGetN( Pptr ),
									mxGetM( gPtr ), mxGetN( gPtr ),
									mxGetM( Cptr ), mxGetN( Cptr ),
									mxGetM( cPtr ), mxGetN( cPtr ),
									mxGetM( zLowPtr ), mxGetN( zLowPtr ),
									mxGetM( zUppPtr ), mxGetN( zUppPtr )
									);
		mexErrMsgTxt( "[qpDUNES] Error: Detected c, but dimensions are inconsistent." );
	}
	if ( ( zLow != 0 ) && ( mxGetM( zLowPtr )*mxGetN( zLowPtr ) != nI*nZ+nX ) )	{
		mexPrintf( "[qpDUNES] Dimensions of data received:\n  nI  : %d\n  H   : %d x %d\n  P   : %d x %d\n  g   : %d x %d\n  C   : %d x %d\n  c   : %d x %d\n  zLow: %d x %d\n  zUpp: %d x %d\n",
									nI,
									mxGetM( Hptr ), mxGetN( Hptr ),
									mxGetM( Pptr ), mxGetN( Pptr ),
									mxGetM( gPtr ), mxGetN( gPtr ),
									mxGetM( Cptr ), mxGetN( Cptr ),
									mxGetM( cPtr ), mxGetN( cPtr ),
									mxGetM( zLowPtr ), mxGetN( zLowPtr ),
									mxGetM( zUppPtr ), mxGetN( zUppPtr )
									);
		mexErrMsgTxt( "[qpDUNES] Error: Detected zLow, but dimensions are inconsistent." );
	}
	if ( ( zUpp != 0 ) && ( mxGetM( zUppPtr )*mxGetN( zUppPtr ) != nI*nZ+nX ) )	{
		mexPrintf( "[qpDUNES] Dimensions of data received:\n  nI  : %d\n  H   : %d x %d\n  P   : %d x %d\n  g   : %d x %d\n  C   : %d x %d\n  c   : %d x %d\n  zLow: %d x %d\n  zUpp: %d x %d\n",
									nI,
									mxGetM( Hptr ), mxGetN( Hptr ),
									mxGetM( Pptr ), mxGetN( Pptr ),
									mxGetM( gPtr ), mxGetN( gPtr ),
									mxGetM( Cptr ), mxGetN( Cptr ),
									mxGetM( cPtr ), mxGetN( cPtr ),
									mxGetM( zLowPtr ), mxGetN( zLowPtr ),
									mxGetM( zUppPtr ), mxGetN( zUppPtr )
									);
		mexErrMsgTxt( "[qpDUNES] Error: Detected zUpp, but dimensions are inconsistent." );
	}
	if ( D_F != 0 ) {
		mexErrMsgTxt( "[qpDUNES] Error: Affine constraints not yet supported." );
	}
	if ( dLow != 0 ) {
		mexErrMsgTxt( "[qpDUNES] Error: Affine constraints not yet supported." );
	}
	if ( dUpp != 0 ) {
		mexErrMsgTxt( "[qpDUNES] Error: Affine constraints not yet supported." );
	}



	/* III) ACTUALLY SET UP QP42 PROBLEM: */
	qpDUNES_setupMatlab( &qpDataGlobal, nI, nX, nU, nD, optionsPtr );

	/* transpose matrices (Fortran to C) */
	/* for C;
	 * H and P should be symmetric
	 */
	real_t* C = 0;
	C = new real_t[nI*nX*nZ];
	for (uint_t kk=0; kk<nI; ++kk) {
		convertFortranToC( &(C[kk*nX*nZ]), &(C_F[kk*nX*nZ]), nX, nZ );
	}
	/* for D */
	real_t* D = 0;
	/*D = new real_t[nDttl];
	convertFortranToC( D, D_F, nDttl, nZ );*/


	/* combine H and P */
	uint_t cIdx, rIdx;
	real_t* H = 0;
	H = new real_t[nI*nZ*nZ+nX*nX];
	/* parse H and P Fortran-style (column-major); note that H and P are by definition self-transposed */
	for ( cIdx=0; cIdx<nI*nZ; ++cIdx ) {
		for ( rIdx=0; rIdx<nZ; ++rIdx ) {
			H[rIdx+cIdx*nZ] = H_in[rIdx+cIdx*nZ];
		}
	}
	for ( cIdx=0; cIdx<nX; ++cIdx ) {
		for ( rIdx=0; rIdx<nX; ++rIdx ) {
			H[nI*nZ*nZ+rIdx+cIdx*nX] = P[rIdx+cIdx*nX];
		}
	}


	/* setup Data */
	if (qpDataGlobal->options.printLevel >= 3) {
		mexPrintf( "Received QP problem of size [nI = %d, nX = %d, nU = %d]\n", nI, nX, nU );
	}


	/* set up a QP problem */
	statusFlag = qpDUNES_init( qpDataGlobal, H, g, C, c, zLow,zUpp, D, dLow, dUpp );
	if ( statusFlag != QPDUNES_OK ) {
		mexPrintf( "Problem setup returned error code %d", statusFlag );
		mexErrMsgTxt( "[qpDUNES] Error: Problem setup failed!" );
	}


	if( nD ) delete[] nD;
	if( H ) delete[] H;
	if( C ) delete[] C;
	if( D ) delete[] D;


	return;
}
/*<<< END OF initQP */



/*
 *	u p d a t e A l l I n t e r v a l s
 */
void updateAllIntervals( qpData_t* qpData, int nrhs, const mxArray* const prhs[] )
{
	/* 0) VARIABLE DECLARATIONS: */
	return_t statusFlag;

	mxArray* Hptr = 0;
	real_t* H_in = 0;
	mxArray* Pptr = 0;
	real_t* P = 0;
	mxArray* gPtr = 0;
	real_t* g = 0;

	mxArray* Cptr = 0;
	real_t* C_F = 0;
	mxArray* cPtr = 0;
	real_t* c = 0;

	mxArray* zLowPtr = 0;
	real_t* zLow = 0;
	mxArray* zUppPtr = 0;
	real_t* zUpp = 0;

	mxArray* Dptr = 0;
	real_t* D_F = 0;
	mxArray* dLowPtr = 0;
	real_t* dLow = 0;
	mxArray* dUppPtr = 0;
	real_t* dUpp = 0;

	const mxArray* optionsPtr = 0;


	/* I) GET OUTPUT POINTERS: */
	Hptr = (mxArray*)prhs[0];	/* to get/check dimensions */
	Pptr = (mxArray*)prhs[1];
	gPtr = (mxArray*)prhs[2];
	Cptr = (mxArray*)prhs[3];
	cPtr = (mxArray*)prhs[4];
	zLowPtr = (mxArray*)prhs[5];
	zUppPtr = (mxArray*)prhs[6];
	Dptr = (mxArray*)prhs[7];
	dLowPtr = (mxArray*)prhs[8];
	dUppPtr = (mxArray*)prhs[9];
	/* exatract data */
	H_in = (real_t*) mxGetPr( Hptr );
	P = (real_t*) mxGetPr( Pptr );
	g = (real_t*) mxGetPr( gPtr );
	C_F = (real_t*) mxGetPr( Cptr );
	c = (real_t*) mxGetPr( cPtr );
	zLow = (real_t*) mxGetPr( zLowPtr );
	zUpp = (real_t*) mxGetPr( zUppPtr );
	D_F = (real_t*) mxGetPr( Dptr );
	dLow = (real_t*) mxGetPr( dLowPtr );
	dUpp = (real_t*) mxGetPr( dUppPtr );




	/* II) CONSISTENCY CHECKS: */
	uint_t nI = qpData->nI;
	uint_t nX = qpData->nX;
	uint_t nZ = qpData->nZ;

	/* 4) check dimension consistency */
	if (qpData->options.printLevel >= 3) {
		mexPrintf( "[qpDUNES] Dimensions of data received:\n  nI  : %d\n  H   : %d x %d\n  P   : %d x %d\n  g   : %d x %d\n  C   : %d x %d\n  c   : %d x %d\n  zLow: %d x %d\n  zUpp: %d x %d\n",
							nI,
							mxGetM( Hptr ), mxGetN( Hptr ),
							mxGetM( Pptr ), mxGetN( Pptr ),
							mxGetM( gPtr ), mxGetN( gPtr ),
							mxGetM( Cptr ), mxGetN( Cptr ),
							mxGetM( cPtr ), mxGetN( cPtr ),
							mxGetM( zLowPtr ), mxGetN( zLowPtr ),
							mxGetM( zUppPtr ), mxGetN( zUppPtr )
							);
	}
	if ( ( H_in != 0 ) != ( P != 0 ) ) {
		mexErrMsgTxt( "[qpDUNES] Error: Either both H and P, or none of them have to be provided for update." );
	}
	if ( ( H_in != 0 ) && ( ( mxGetM( Hptr ) != nZ ) || ( mxGetN( Hptr ) != nI*nZ ) ) )	{
		mexErrMsgTxt( "[qpDUNES] Error: Detected H, but dimensions are inconsistent." );
	}
	if ( ( P != 0 ) && ( ( mxGetM( Pptr ) != nX ) || ( mxGetN( Pptr ) != nX ) ) ) {
		mexErrMsgTxt( "[qpDUNES] Error: Detected P, but dimensions are inconsistent." );
	}
	if ( ( g != 0 ) && ( mxGetM( gPtr )*mxGetN( gPtr ) != nI*nZ + nX ) ) {
		mexErrMsgTxt( "[qpDUNES] Error: Detected g, but dimensions are inconsistent." );
	}
	if ( ( C_F != 0 ) && ( ( mxGetM( Cptr ) != nX ) || ( mxGetN( Cptr ) != nI*nZ ) ) ) {
		mexErrMsgTxt( "[qpDUNES] Error: Detected C, but dimensions are inconsistent." );
	}
	if ( ( c != 0 ) && ( mxGetM( cPtr )*mxGetN( cPtr ) != nI*nX ) ) {
		mexErrMsgTxt( "[qpDUNES] Error: Detected c, but dimensions are inconsistent." );
	}
	if ( ( zLow != 0 ) && ( mxGetM( zLowPtr )*mxGetN( zLowPtr ) != nI*nZ+nX ) )	{
		mexErrMsgTxt( "[qpDUNES] Error: Detected zLow, but dimensions are inconsistent." );
	}
	if ( ( zUpp != 0 ) && ( mxGetM( zUppPtr )*mxGetN( zUppPtr ) != nI*nZ+nX ) )	{
		mexErrMsgTxt( "[qpDUNES] Error: Detected zUpp, but dimensions are inconsistent." );
	}
	if ( D_F != 0 ) {
		mexErrMsgTxt( "[qpDUNES] Error: Affine constraints not yet supported." );
	}
	if ( dLow != 0 ) {
		mexErrMsgTxt( "[qpDUNES] Error: Affine constraints not yet supported." );
	}
	if ( dUpp != 0 ) {
		mexErrMsgTxt( "[qpDUNES] Error: Affine constraints not yet supported." );
	}



	/* III) ACTUALLY SET UP QP42 PROBLEM: */
	/* transpose matrices (Fortran to C) */
	/* only for C, D
	 * H should be symmetric
	*/
	real_t* C = 0;
	if (C_F != 0) {
		C = new real_t[nI*nX*nZ];
		for (uint_t kk=0; kk<nI; ++kk) {
			convertFortranToC( &(C[kk*nX*nZ]), &(C_F[kk*nX*nZ]), nX, nZ );
		}
	}
	real_t* D = 0;
	if (D_F != 0) {
		/* TODO: transpose D */
	}

	/* combine H and P (only given if LTV system) */
	uint_t cIdx, rIdx;
	real_t* H = 0;
	if( (H_in != 0) && (P != 0) ) {
		H = new real_t[nI*nZ*nZ+nX*nX];
		/* parse H and P Fortran-style (column-major); note that H and P are by definition self-transposed */
		for ( cIdx=0; cIdx<nI*nZ; ++cIdx ) {
			for ( rIdx=0; rIdx<nZ; ++rIdx ) {
				H[rIdx+cIdx*nZ] = H_in[rIdx+cIdx*nZ];
			}
		}
		for ( cIdx=0; cIdx<nX; ++cIdx ) {
			for ( rIdx=0; rIdx<nX; ++rIdx ) {
				H[nI*nZ*nZ+rIdx+cIdx*nX] = P[rIdx+cIdx*nX];
			}
		}
	}


	statusFlag = qpDUNES_updateData( qpData, H, g, C, c, zLow, zUpp, D, dLow, dUpp );	/* TODO: support general constraints	*/
	if (statusFlag != QPDUNES_OK) {
		mexErrMsgTxt( "QP data update failed." );
		return;
	}


	/* print which data fields are updated */
	if (qpData->options.printLevel >= 3) {
		mexPrintf( "Received Data updates for QP of size [nI = %d, nX = %d, nU = %d]:\n", nI, nX, qpData->nU );
		( H_in != 0 ) ? mexPrintf( "H\n" ) : mexPrintf( "" );
		( P != 0 ) ? mexPrintf( "P\n" ) : mexPrintf( "" );
		( g != 0 ) ? mexPrintf( "g\n" ) : mexPrintf( "" );
		( C_F != 0 ) ? mexPrintf( "C\n" ) : mexPrintf( "" );
		( c != 0 ) ? mexPrintf( "c\n" ) : mexPrintf( "" );
		( zLow != 0 ) ? mexPrintf( "zLow\n" ) : mexPrintf( "" );
		( zUpp != 0 ) ? mexPrintf( "zUpp\n" ) : mexPrintf( "" );
		( D_F != 0 ) ? mexPrintf( "D\n" ) : mexPrintf( "" );
		( dLow != 0 ) ? mexPrintf( "dLow\n" ) : mexPrintf( "" );
		( dUpp != 0 ) ? mexPrintf( "dUpp\n" ) : mexPrintf( "" );
	}


	if( H ) delete[] H;
	if( C ) delete[] C;
	if( D ) delete[] D;


	return;
}
/*<<< END OF updateAllIntervals */


/*
 *	u p d a t e S i n g l e I n t e r v a l
 */
void updateSingleInterval( qpData_t* qpData, int nrhs, const mxArray* const prhs[] )
{
	/* 0) VARIABLE DECLARATIONS: */
	return_t statusFlag;

	mxArray* IidxPtr = 0;
	real_t* IidxFloat = 0;

	mxArray* Hptr = 0;
	real_t* H = 0;
	mxArray* gPtr = 0;
	real_t* g = 0;

	mxArray* Cptr = 0;
	real_t* C_F = 0;
	mxArray* cPtr = 0;
	real_t* c = 0;

	mxArray* zLowPtr = 0;
	real_t* zLow = 0;
	mxArray* zUppPtr = 0;
	real_t* zUpp = 0;

	mxArray* Dptr = 0;
	real_t* D_F = 0;
	mxArray* dLowPtr = 0;
	real_t* dLow = 0;
	mxArray* dUppPtr = 0;
	real_t* dUpp = 0;

	const mxArray* optionsPtr = 0;


	/* I) GET OUTPUT POINTERS: */
	IidxPtr = (mxArray*)prhs[0];	/* to get/check dimensions */
	Hptr = (mxArray*)prhs[1];
	gPtr = (mxArray*)prhs[2];
	Cptr = (mxArray*)prhs[3];
	cPtr = (mxArray*)prhs[4];
	zLowPtr = (mxArray*)prhs[5];
	zUppPtr = (mxArray*)prhs[6];
	Dptr = (mxArray*)prhs[7];
	dLowPtr = (mxArray*)prhs[8];
	dUppPtr = (mxArray*)prhs[9];
	/* exatract data */
	IidxFloat = (real_t*) mxGetPr( IidxPtr );
	H = (real_t*) mxGetPr( Hptr );
	g = (real_t*) mxGetPr( gPtr );
	C_F = (real_t*) mxGetPr( Cptr );
	c = (real_t*) mxGetPr( cPtr );
	zLow = (real_t*) mxGetPr( zLowPtr );
	zUpp = (real_t*) mxGetPr( zUppPtr );
	D_F = (real_t*) mxGetPr( Dptr );
	dLow = (real_t*) mxGetPr( dLowPtr );
	dUpp = (real_t*) mxGetPr( dUppPtr );




	/* II) CONSISTENCY CHECKS: */
	uint_t nI = qpData->nI;
	uint_t nX = qpData->nX;
	uint_t nZ = qpData->nZ;

	/* get index of interval to be updated */
	if ( mxGetM( IidxPtr )*mxGetN( IidxPtr ) != 1 )	 mexErrMsgTxt( "[qpDUNES] Error: Index of interval to be updated needs to be specified." );
	int_t Iidx = (int_t)(*IidxFloat);
	if ( (Iidx < 0) || (Iidx > nI) )  mexErrMsgTxt( "[qpDUNES] Error: Interval index needs to be between 0 and nI." );
	uint_t nV = qpData->intervals[Iidx]->nV;


	/* check dimension consistency */
	if (qpData->options.printLevel >= 3) {
		mexPrintf( "[qpDUNES] Dimensions received:\n  Iidx  : %d\n  H   : %d x %d\n  g   : %d x %d\n  C   : %d x %d\n  c   : %d x %d\n  zLow: %d x %d\n  zUpp: %d x %d\n",
							nI,
							mxGetM( Hptr ), mxGetN( Hptr ),
							mxGetM( gPtr ), mxGetN( gPtr ),
							mxGetM( Cptr ), mxGetN( Cptr ),
							mxGetM( cPtr ), mxGetN( cPtr ),
							mxGetM( zLowPtr ), mxGetN( zLowPtr ),
							mxGetM( zUppPtr ), mxGetN( zUppPtr )
							);
	}
	if ( ( H != 0 ) && ( ( mxGetM( Hptr ) != nV ) || ( mxGetN( Hptr ) != nV ) ) )	{
		mexErrMsgTxt( "[qpDUNES] Error: Detected H, but dimensions are inconsistent." );
	}
	if ( ( g != 0 ) && ( mxGetM( gPtr )*mxGetN( gPtr ) != nV ) ) {
		mexErrMsgTxt( "[qpDUNES] Error: Detected g, but dimensions are inconsistent." );
	}

	if ( ( Iidx == nI ) && ( C_F != 0 )) {
		mexErrMsgTxt( "[qpDUNES] Error: Terminal interval specified, but C detected." );
	}
	if ( ( Iidx == nI ) && ( c != 0 )) {
		mexErrMsgTxt( "[qpDUNES] Error: Terminal interval specified, but c detected." );
	}
	if ( ( C_F != 0 ) && ( ( mxGetM( Cptr ) != nX ) || ( mxGetN( Cptr ) != nZ ) ) ) {
		mexErrMsgTxt( "[qpDUNES] Error: Detected C, but dimensions are inconsistent." );
	}
	if ( ( c != 0 ) && ( mxGetM( cPtr )*mxGetN( cPtr ) != nV ) ) {
		mexErrMsgTxt( "[qpDUNES] Error: Detected c, but dimensions are inconsistent." );
	}

	if ( ( zLow != 0 ) && ( mxGetM( zLowPtr )*mxGetN( zLowPtr ) != nV ) )	{
		mexErrMsgTxt( "[qpDUNES] Error: Detected zLow, but dimensions are inconsistent." );
	}
	if ( ( zUpp != 0 ) && ( mxGetM( zUppPtr )*mxGetN( zUppPtr ) != nV ) )	{
		mexErrMsgTxt( "[qpDUNES] Error: Detected zUpp, but dimensions are inconsistent." );
	}

	if ( D_F != 0 ) {
		mexErrMsgTxt( "[qpDUNES] Error: Affine constraints not yet supported." );
	}
	if ( dLow != 0 ) {
		mexErrMsgTxt( "[qpDUNES] Error: Affine constraints not yet supported." );
	}
	if ( dUpp != 0 ) {
		mexErrMsgTxt( "[qpDUNES] Error: Affine constraints not yet supported." );
	}



	/* III) ACTUALLY SET UP QP42 PROBLEM: */
	/* transpose matrices (Fortran to C) */
	/* only for C, D
	 * H should be symmetric
	*/
	real_t* C = 0;
	if (C_F != 0) {
		C = new real_t[nI*nX*nZ];
		for (uint_t kk=0; kk<nI; ++kk) {
			convertFortranToC( &(C[kk*nX*nZ]), &(C_F[kk*nX*nZ]), nX, nZ );
		}
	}
	real_t* D = 0;
	if (D_F != 0) {
		/* TODO: transpose D */
	}


	statusFlag = qpDUNES_updateIntervalData( qpData, qpData->intervals[Iidx], H, g, C, c, zLow,zUpp, D,dLow,dUpp, 0 );
	if (statusFlag != QPDUNES_OK) {
		mexErrMsgTxt( "QP data update failed." );
		return;
	}

	if( C ) delete[] C;
	if( D ) delete[] D;

	return;
}
/*<<< END OF updateSingleInterval */



/*
 *	s o l v e M a t l a b
 */
void solveMatlab( int nlhs, mxArray* plhs[], int nrhs, const mxArray* const prhs[] )
{
	/* GET INPUTS */
//	mxArray* x0Ptr = (mxArray*)prhs[0];
//	real_t* x0 = (real_t*) mxGetPr( x0Ptr );


	/* CONSISTENCY CHECKS */
	/* 1)   Check for proper number of output arguments. */
	if (qpDataGlobal->options.logLevel == QPDUNES_LOG_ALL_DATA ) {
		if ( nlhs > 6 )	mexErrMsgTxt( "[qpDUNES] Error: At most six output arguments are allowed: \n                 [zOpt, status, lambda, mu, objVal, log]!" );
	}
	else {
		if ( nlhs > 5 )	mexErrMsgTxt( "[qpDUNES] Error: At most five output arguments are allowed: \n                 [zOpt, status, lambda, mu, objVal]!" );
	}
	if ( nlhs < 1 )	mexErrMsgTxt( "[qpDUNES] Error: At least one output argument is required: [zOpt,...]!" );


	/* ALLOCATE OUTPUTS */
	allocateOutputsQP( plhs,nlhs, qpDataGlobal->nI,qpDataGlobal->nX,qpDataGlobal->nZ );


	/* SOLVE QP42 PROBLEM: */
	return_t statusFlag = qpDUNES_solve( qpDataGlobal );
	if ( statusFlag != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND ) {
		mexPrintf( "qpDUNES returned flag %d\n", statusFlag );
		if (qpDataGlobal->options.logLevel == QPDUNES_LOG_ALL_DATA ) {
			if ( nlhs == 5 ) {
				fullLogging( qpDataGlobal, &(plhs[4]) );
			}
		}
		mexErrMsgTxt( "[qpDUNES] Error: Problem could not be solved!" );
	}


	/* V) PASS SOLUTION ON TO MATLAB: */
	obtainOutputsQP( qpDataGlobal, plhs, nlhs, statusFlag );


	/* VI) PASS DETAILED LOG INFORMATION ON TO MATLAB: */
	if (qpDataGlobal->options.logLevel == QPDUNES_LOG_ALL_DATA ) {
		if ( nlhs == 5 ) {
			fullLogging( qpDataGlobal, &(plhs[4]) );
		}
	}


	return;
}
/*<<< END OF solveMatlab */




/*
 *	m e x F u n c t i o n
 */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
	/* I) CONSISTENCY CHECKS: Ensure that first input is a string */
	if ( mxIsChar( prhs[0] ) != 1 )
		mexErrMsgTxt( "[qpDUNES] Error: First input argument must be a string!" );

	char* typeString = mxArrayToString( prhs[0] );


	/* SELECT QPDUNES ACTION TO BE PERFORMED */
	/* a) initial setup */
	if ( ( strcmp( typeString,"i" ) == 0 ) ||
		 ( strcmp( typeString,"I" ) == 0 ) ||
		 ( strcmp( typeString,"init" ) == 0 ) ||
		 ( strcmp( typeString,"Init" ) == 0 ) ||
		 ( strcmp( typeString,"INIT" ) == 0 ) )
	{
		if (nlhs > 0) {
			mexErrMsgTxt( "[qpDUNES] Error: Too many output arguments for action 'init'.\n                 Type 'help qpDUNES' for further information.");
		}

		switch ( nrhs-1 ) {
			case 11:
			case 12:
				initQP( nrhs-1, prhs+1 );	/* shift input pointers to begin of data */
				break;

			default:
				mexPrintf( "Received %d input arguments.\n", nrhs );
				mexErrMsgTxt( "[qpDUNES] Error: Wrong number of input arguments for action 'init'.\n                 Type 'help qpDUNES' for further information." );
				break;
		}

		return;
	}


	/* b) full data update */
	if ( ( strcmp( typeString,"u" ) == 0 ) ||
		 ( strcmp( typeString,"U" ) == 0 ) ||
		 ( strcmp( typeString,"update" ) == 0 ) ||
		 ( strcmp( typeString,"Update" ) == 0 ) ||
		 ( strcmp( typeString,"UPDATE" ) == 0 ) )
	{
		if( qpDataGlobal == 0 )
			mexErrMsgTxt( "[qpDUNES] Error: QP data needs to initialized at least once before action 'update' can be performed.\n                 Type 'help qpDUNES' for further information." );

		if (nlhs > 0) {
			mexErrMsgTxt( "[qpDUNES] Error: Too many output arguments for action 'update'.\n                 Type 'help qpDUNES' for further information.");
		}

		switch ( nrhs-1 ) {
			case 10:
				updateAllIntervals( qpDataGlobal, nrhs-1, prhs+1 );	/* shift input pointers to begin of data */
				break;

			default:
				mexErrMsgTxt( "[qpDUNES] Error: Wrong number of input arguments for action 'update'.\n                 Type 'help qpDUNES' for further information.");
		}

		return;
	}


	/* c) single interval update */
	if ( ( strcmp( typeString,"ui" ) == 0 ) ||
		 ( strcmp( typeString,"iu" ) == 0 ) ||
		 ( strcmp( typeString,"UI" ) == 0 ) ||
		 ( strcmp( typeString,"IU" ) == 0 ) ||
		 ( strcmp( typeString,"updateinterval" ) == 0 ) ||
		 ( strcmp( typeString,"updateInterval" ) == 0 ) ||
		 ( strcmp( typeString,"Updateinterval" ) == 0 ) ||
		 ( strcmp( typeString,"UpdateInterval" ) == 0 ) ||
		 ( strcmp( typeString,"UPDATEINTERVAL" ) == 0 ) ||
		 ( strcmp( typeString,"intervalupdate" ) == 0 ) ||
		 ( strcmp( typeString,"intervalUpdate" ) == 0 ) ||
		 ( strcmp( typeString,"Intervalupdate" ) == 0 ) ||
		 ( strcmp( typeString,"IntervalUpdate" ) == 0 ) ||
		 ( strcmp( typeString,"INTERVALUPDATE" ) == 0 ) )
	{
		if( qpDataGlobal == 0 )
			mexErrMsgTxt( "[qpDUNES] Error: QP data needs to initialized at least once before action 'updateInterval' can be performed.\n                 Type 'help qpDUNES' for further information." );

		if (nlhs > 0) {
			mexErrMsgTxt( "[qpDUNES] Error: Too many output arguments for action 'updateInterval'.\n                 Type 'help qpDUNES' for further information.");
		}

		switch ( nrhs-1 ) {
			case 10:
				updateSingleInterval( qpDataGlobal, nrhs-1, prhs+1 );	/* shift input pointers to begin of data */
				break;

			default:
				mexErrMsgTxt( "[qpDUNES] Error: Wrong number of input arguments for action 'updateInterval'.\n                 Type 'help qpDUNES' for further information.");
		}

		return;
	}


	/* d) solve */
	if ( ( strcmp( typeString,"s" ) == 0 ) ||
		 ( strcmp( typeString,"S" ) == 0 ) ||
		 ( strcmp( typeString,"solve" ) == 0 ) ||
		 ( strcmp( typeString,"Solve" ) == 0 ) ||
		 ( strcmp( typeString,"SOLVE" ) == 0 ) )
	{
		if( qpDataGlobal == 0 )
			mexErrMsgTxt( "[qpDUNES] Error: QP data needs to initialized before action 'solve' can be performed.\n                 Type 'help qpDUNES' for further information." );

		if (nrhs > 2) {
			mexErrMsgTxt( "[qpDUNES] Error: Too many input arguments for action 'solve'.\n                 Type 'help qpDUNES' for further information.");
		}

		solveMatlab( nlhs, plhs, nrhs-1, prhs+1 );

		return;
	}


	/* e) shift */
	if ( ( strcmp( typeString,"sh" ) == 0 ) ||
		 ( strcmp( typeString,"SH" ) == 0 ) ||
		 ( strcmp( typeString,"Sh" ) == 0 ) ||
		 ( strcmp( typeString,"shift" ) == 0 ) ||
		 ( strcmp( typeString,"Shift" ) == 0 ) ||
		 ( strcmp( typeString,"SHIFT" ) == 0 ) )
	{
		if( qpDataGlobal == 0 )
			mexErrMsgTxt( "[qpDUNES] Error: QP data needs to initialized before action 'shift' can be performed.\n                 Type 'help qpDUNES' for further information." );

		if (nlhs > 0) {
			mexErrMsgTxt( "[qpDUNES] Error: Too many output arguments for action 'shift'.\n                 Type 'help qpDUNES' for further information.");
		}

		if (nrhs-1 > 0) {
			mexErrMsgTxt( "[qpDUNES] Error: Too many input arguments for action 'shift'.\n                 Type 'help qpDUNES' for further information.");
		}


		/* shift to prepare QP for next solution */
		qpDUNES_shiftLambda( qpDataGlobal );			/* shift multipliers */
		qpDUNES_shiftIntervals( qpDataGlobal );			/* shift intervals (particulary important when using qpOASES for underlying local QPs) */


		return;
	}


	mexErrMsgTxt( "[qpDUNES] Error: Undefined first input argument!\n                 Type 'help qpDUNES' for further information." );

	return;
}

/*
 *	end of file
 */
