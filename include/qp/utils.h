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
 *	License along with qp42; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qp/utils.h
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */


#ifndef QPDUNES_UTILS_H
#define QPDUNES_UTILS_H


#include <stdio.h>
#include <stdarg.h>

#ifndef WIN32
    #include <sys/time.h>
#endif

#include <qp/types.h>

#ifdef __QPDUNES_PARALLEL__
	#include <omp.h>
#endif

#ifdef __MATLAB__
	#include "mex.h"
#endif


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_free(	real_t** data
				);

/** 
 *	\brief A bit safer calloc. Internally does assert against NULL pointer.
 *
 *	\author Milan Vukov
 *	\version 1.0beta
 *	\date 2014
 */
void* qpDUNES_calloc(size_t num, size_t size);


/**
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_intFree(	int_t** data
					);



/**
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
/*inline real_t* offsetArray(	real_t* data,*/
const real_t* offsetArray(	const real_t* const data,
							int_t offset
							);




/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
/*return_t qp42_setupMatrix(	matrix_t* const to,
							const real_t* const from,
							int_t nRows,
							int_t nCols
							);*/
sparsityType_t qpDUNES_detectMatrixSparsity(	const real_t* const M,
											int_t nRows,
											int_t nCols
											);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_updateMatrixData(	matrix_t* const to,
								const real_t* const from,
								int_t nRows,
								int_t nCols
								);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_setupZeroMatrix(	int_t nRows,
								int_t nCols,
								matrix_t* to
								);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_setMatrixNull(	matrix_t* const matrix
								);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_existsMatrix(	matrix_t* matrix
							);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_existsVector(	vector_t* vector
							);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_setupIdentityMatrix(	matrix_t* to
									);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_setupScaledIdentityMatrix(	int_t nRows,
											real_t scalar,
											matrix_t* to
											);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_setupVector(	vector_t* const to,
							const real_t* const from,
							int_t n
							);

/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_updateVector(	vector_t* const to,
							const real_t* const from,
							int_t n
							);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_updateSimpleBoundVector(	qpData_t* qpData,
										vector_t* const to,
										const real_t* const dBnd,
										const real_t* const xBnd,
										const real_t* const uBnd
										);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_updateConstraintVector( 	vector_t* const to,
										const real_t* const dBnd,
										int_t nD
										);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_setupZeroVector(	vector_t* const to,
								int_t n
								);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_setupUniformVector(	vector_t* const to,
									real_t value,
									int_t n
									);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_copyVector(	vector_t* const to,
							const vector_t* const from,
							int_t n
							);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
/*return_t qp42_copyIntVector(	intVector_t* const to,
								const intVector_t* const from,
								int_t n
								);*/


/**
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_copyMatrix(	matrix_t* const to,
							const matrix_t* const from,
							int_t dim0,
							int_t dim1
							);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
return_t qpDUNES_makeMatrixDense( matrix_t* const M, 
								int_t dim0,
								int_t dim1
								);


return_t qpDUNES_transposeMatrix(	matrix_t* const to,
								const matrix_t* const from,
								int_t dim0,
								int_t dim1
								);


return_t qpDUNES_selftransposeMatrix(	matrix_t* const M,
									int_t dim			/**< leading and secondary dimension of M */
									);


return_t qpDUNES_copyArray(	real_t* const to,
							const real_t* const from,
							int_t n
							);



/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
/* return_t getLambdaOffset(	xn_vector_t* lambda,
 							int_t idx
 							);*/




/** max routine */
/*inline int_t qp42_max(	int_t a, */
int_t qpDUNES_max(	int_t a,
						int_t b );


/** min routine */
/*inline int_t qp42_min(	int_t a, */
int_t qpDUNES_min(	int_t a,
						int_t b );


/** real_t max routine */
/*inline real_t qp42_fmax(	real_t a, */
real_t qpDUNES_fmax(	real_t a,
							real_t b );


/** real_t min routine  */
/*inline real_t qp42_fmin(	real_t a, */
real_t qpDUNES_fmin(	real_t a,
							real_t b );


/** sign routine  */
int_t qpDUNES_sign(	const qpData_t* const qpData,
						real_t a
						);



/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
/* inline void qpDUNES_assertOK(	return_t statusFlag, */
void qpDUNES_assertOK(	return_t statusFlag,
							char* fileName,
							int_t lineNumber,
							char* errString 
);




/** 
 *	\brief Prepare results struct
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
/*void qp42_getResults( 	qpData_t qpData,
						qpResults_t results
						);*/



/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
real_t getTime( );



/** 
 *	\brief Customizable low-level printing routine
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printStrArgs(	const char* const string,
						... 
						);


/** 
 *	\brief Low-level print to file
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printStrArgsToFile(	FILE* filePtr,
								const char* const string,
								...
								);


/**
 *	\brief Customizable low-level printing routine
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printStrArgsList(	const char* const string,
							va_list printArgs 
							);


/** 
 *	\brief Low-level print to file
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printStrArgsListToFile(	FILE* filePtr,
									const char* const string,
									va_list printArgs
									);


/**
 *	\brief Customizable printf routine
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printf(	const char* const string,
					... 
					);





/** 
 *	\brief Customizable printf routine
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printf_noNewLine(	const char* const string,
							...
							);


/**
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printSuccess( const qpData_t* const qpData,
						const char* const string,
						...
						);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printWarning(	const qpData_t* const qpData,
						const char* const fileName,
						const int_t lineNumber,
						const char* const errString );


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printError(	const qpData_t* const qpData,
						const char* const fileName,
						const int_t lineNumber,
						const char* const errString,
						...
						);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printDebugInfo( const char* const string );


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printMatrixData(	const real_t* const M,
							const int_t dim0,
							const int_t dim1,
							const char* const string,
							...
							);


/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printMatrixDataToFile(	const real_t* const M,
									const int_t dim0,
									const int_t dim1,
									char* fileName,
									const char* const varNameString,
									...
									);



/**
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printNewtonHessian(	const qpData_t* const qpData,
								const xn2x_matrix_t* const hessian
								);


/**
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printNewtonHessianToFile(	const qpData_t* const qpData,
									const xn2x_matrix_t* const hessian,
									const char* const fileName,
									const char* const variableName,
									...
									);


/**
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
void qpDUNES_printCholNewtonHessian(	const qpData_t* const qpData,
									const xn2x_matrix_t* const hessian
									);


/**
 *	\brief qpDUNES header information
 *
 *	\author Janick Frasch
 *	\version 1.0beta
 *	\date 2013
 */
void qpDUNES_printHeader( qpData_t* qpData );

/**
 *	\brief Get error string
 *
 *	\author Milan Vukov
 *	\version 1.0beta
 *	\date 2014
 */
const char* qpDUNES_getErrorString( int code );


#endif	/* QPDUNES_UTILS_H */


/*
 *	end of file
 */
