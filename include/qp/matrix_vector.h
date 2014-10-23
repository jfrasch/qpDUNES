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
 *	\file include/qp/matrix_vector.h
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */


#ifndef QP42_MATRIX_VECTOR_H
#define QP42_MATRIX_VECTOR_H


#include <qp/types.h>
#include <qp/qpdunes_utils.h>
#include <math.h>


return_t multiplyQx(	qpData_t* const qpData,
						x_vector_t* const res,
						const xx_matrix_t* const Q,
						const x_vector_t* const x 	);


return_t multiplyInvQx(	qpData_t* const qpData,
						x_vector_t* const res,
						const xx_matrix_t* const cholQ,
						const x_vector_t* const y 	);


return_t multiplyRu(	qpData_t* const qpData,
						u_vector_t* const res,
						const uu_matrix_t* const R,
						const u_vector_t* const u 	);


return_t multiplyInvRu(	qpData_t* const qpData,
						u_vector_t* const res,
						const uu_matrix_t* const cholR,
						const u_vector_t* const y 	);


real_t multiplyzHz(	qpData_t* const qpData,
					const vv_matrix_t* const H,
					const z_vector_t* const z,
					 const int_t nV 	);


return_t multiplyInvHz(	qpData_t* const qpData,
						z_vector_t* const res,
						const vv_matrix_t* const cholH,
						const z_vector_t* const z,
						const int_t nV 	);


return_t multiplyCz(	qpData_t* const qpData,
						x_vector_t* const res,
						const xz_matrix_t* const C,
						const z_vector_t* const z 	);




return_t multiplyCTy(	qpData_t* const qpData,
						z_vector_t* const res,
						const xz_matrix_t* const C,
						const x_vector_t* const y 	);


/** Inverse matrix times matrix product res = Q^-1 * A */
return_t multiplyAInvQ(	qpData_t* const qpData,
						xx_matrix_t* const res,
						const xx_matrix_t* const C,
						const vv_matrix_t* const cholH
						);

return_t multiplyInvQAT(	qpData_t* const qpData,
						xx_matrix_t* const res,
						const xx_matrix_t* const cholQ,
						const xx_matrix_t* const A,
						x_vector_t* const vecTmp
						);


return_t getInvQ(	qpData_t* const qpData,
					xx_matrix_t* const res,
					const xx_matrix_t* const cholM1,
					int_t nV
					);


return_t factorizeH( 	qpData_t* const qpData,
						vv_matrix_t* const cholH,
						const vv_matrix_t* const H,
						int_t nV
						);


/* compute x_k+1 = A*x_k + B*u_k * c*/
return_t computeXk1(	qpData_t* const qpData,
						x_vector_t* const res,
						const xx_matrix_t* const A,
						const xu_matrix_t* const B,
						const x_vector_t* const c,
						const x_vector_t* const x,
						const u_vector_t* const u 	);


/* todo: check for consistency with .c file!! */
return_t addScaledLambdaStep(	qpData_t* const qpData,
								xn_vector_t* const res,
								real_t scalingFactor,
								const xn_vector_t* const deltaLambda	);


return_t copyScaleVector(	vector_t* const res,
							real_t scalingFactor,
							const vector_t* const vec,
							int_t len
							);


return_t scaleVector(	vector_t* const res,
						real_t scalingFactor,
						int_t len
						);


return_t addScaledVector(	xn_vector_t* const res,
							real_t scalingFactor,
							const xn_vector_t* const update,
							int_t len
							);


return_t addVectorScaledVector(	vector_t* const res,
								const vector_t* const x,
								real_t scalingFactor,
								const vector_t* const y,
								int_t len
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
return_t qp42_newtonGradientNorm(	qpData_t* qpData,
									xn_vector_t* vec
);


return_t backsolveDenseL(	qpData_t* const qpData,
							real_t* const res, 
							const real_t* const L, 
							const real_t* const b, 
							boolean_t transposed,
							int_t n
							);


return_t backsolveDiagonal(	qpData_t* const qpData,
							real_t* const res, 
							const real_t* const M, 
							const real_t* const b, 
							int_t n
							);


return_t backsolveMatrixDenseDenseL( qpData_t* const qpData,
									real_t* const res,
									 const real_t* const L,
									 const real_t* const M,
									 real_t* const sums,			/**< memory for saving intermediate results (for speedup) */
									boolean_t transposed,
									int_t dim0,
									int_t dim1
									);


return_t backsolveMatrixDenseDenseTL( qpData_t* const qpData,
									 real_t* const res,
									 const real_t* const L,
									 const real_t* const M,
									 real_t* const sums,			/**< memory for saving intermediate results (for speedup) */
									boolean_t transposed,
									int_t dim0,			/* leading dimension of M */
									int_t dim1 			/* secondary dimension of M */
									);


return_t backsolveMatrixDenseIdentityL( qpData_t* const qpData,
										real_t* const res,
										const real_t* const L,
										real_t* const sums,			/**< memory for saving intermediate results (for speedup) */
										int_t dim0 );


return_t backsolveMatrixDiagonalDense( qpData_t* const qpData,
									   real_t* const res,
									   const real_t* const M1,
									   const real_t* const M2,
									   int_t dim0,
									   int_t dim1 );

return_t backsolveMatrixDiagonalDenseT( qpData_t* const qpData,
									   real_t* const res,
									   const real_t* const L,
									   const real_t* const M,
									int_t dim0
									);


return_t backsolveMatrixDiagonalDiagonal( qpData_t* const qpData,
										  real_t* const res,
										  const real_t* const M1,
										  const real_t* const M2,
										 int_t dim0 
										 );


return_t backsolveMatrixDiagonalIdentity( qpData_t* const qpData,
										  real_t* const res,
										  const real_t* const M1,
										  int_t dim0 );


/* ----------------------------------------------
 * Matrix backsolve for L, M both dense
 * compute res for L*res = M^T
 * 
 >>>>>                                            */
return_t backsolveMatrixTDenseDenseL( qpData_t* const qpData,
									  real_t* const res,
									  const real_t* const L,
									  const real_t* const M,			/**< untransposed M */
									  real_t* const sums,				/**< memory for saving intermediate results (for speedup) */
									  boolean_t transposedL,
									  int_t dim0,						/**< dimensions of M */
									  int_t dim1 );


return_t backsolveMatrixTDiagonalDense( qpData_t* const qpData,
										real_t* const res,
										const real_t* const M1,
										const real_t* const M2,	/**< untransposed M2 */
										int_t dim0,				/**< dimensions of M2 */
										int_t dim1 );


return_t multiplyMatrixVector( vector_t* const res,
								const matrix_t* const M,
								const vector_t* const x,
								int_t dim0,
								int_t dim1 	);


real_t multiplyVectorMatrixVector(	const matrix_t* const M,
							   	   	const vector_t* const x,
							   	   	int_t dim0	);


return_t multiplyBlockDiagMatrixVector( vector_t* const res,
										const matrix_t* const M1,
										const matrix_t* const M2,
										const vector_t* const x,
										int_t dimM1,			/**< dimensions of M1 */
										int_t dimM2 			/**< dimensions of M2 */
										);


real_t multiplyBlockDiagVectorMatrixVector( const matrix_t* const M1,
											const matrix_t* const M2,
											const vector_t* const x,
											int_t dimM1,			/**< dimensions of M1 */
											int_t dimM2 			/**< dimensions of M2 */
											);


return_t multiplyMatrixTVector(	vector_t* const res,
								const matrix_t* const M,
								const vector_t* const x,
								int_t dim0,
								int_t dim1 	);


return_t multiplyInvMatrixVector(	qpData_t* const qpData,
									vector_t* const res,
									const matrix_t* const cholH,
									const vector_t* const x,
									int_t dim0						/**< dimension of symmetric matrix */
									);


return_t multiplyInvBlockDiagMatrixVector(	qpData_t* const qpData,
											vector_t* const res,
											const matrix_t* const cholM1,
											const matrix_t* const cholM2,
											const vector_t* const x,
											int_t dimM1,					/**< dimensions of M1 */
											int_t dimM2 					/**< dimensions of M2 */
											);


return_t multiplyInvMatrixMatrix(	qpData_t* const qpData,
									matrix_t* const res,
								  const matrix_t* const cholM1,
								  const matrix_t* const M2,
								  vector_t* const vecTmp,
								  int_t dim0,				/**< leading dimension of A == secondary dimension of A == leading dimension of M2 */
								  int_t dim1				/**< secondary dimension of M2 */
									);


return_t multiplyInvMatrixMatrixT(	qpData_t* const qpData,
									matrix_t* const res,
								  const matrix_t* const cholM1,
								  matrix_t* const M2,
								  vector_t* const vecTmp,
								  int_t dim0,				/**< leading dimension of A == secondary dimension of A == leading dimension of M2 */
								  int_t dim1				/**< secondary dimension of M2 */
									);


return_t multiplyMatrixVectorDense(	real_t* const res,
									const real_t* const M,
									const real_t* const x,
									int_t dim0,
									int_t dim1 	);


real_t multiplyVectorMatrixVectorDense(	const real_t* const M,
										const real_t* const x,
										int_t dim0 );


return_t multiplyMatrixVectorSparse(	real_t* const res,
										const real_t* const M,
										const real_t* const x,
										int_t dim0,
										int_t dim1 	);


return_t multiplyMatrixVectorDiagonal(	real_t* const res,
										const real_t* const M,
										const real_t* const x,
										int_t dim0	);


real_t multiplyVectorMatrixVectorDiagonal(	const real_t* const M,
												const real_t* const x,
												int_t dim0	);


/** Dense generic transposed matrix-vector product res = M.T*x */
return_t multiplyMatrixTVectorDense(	real_t* const res,
										const real_t* const M,	/**< untransposed matrix */
										const real_t* const x,
										int_t dim0,
										int_t dim1		);


/** Sparse generic transposed matrix-vector product res = M.T*x */
return_t multiplyMatrixTVectorSparse(	real_t* const res,
										const real_t* const M,
										const real_t* const x,
										int_t dim0,
										int_t dim1		);


return_t multiplyInvMatrixMatrixDense(	matrix_t* const res,
									const matrix_t* const M1,
									const matrix_t* const M2,
									int_t dim0,						/**< leading dimension of M1 */
									int_t dim1,						/**< secondary dimension of M1 == leading dimension of M2 */
									int_t dim2						/**< secondary dimension of M2 */
									);


/* ----------------------------------------------
 * Dense generic transposed matrix-matrix product
 * res = M1.T*M2
 *
 >>>>>                                            */
void multiplyMatrixTMatrixDenseDense( real_t* const res,
							const real_t* const M1,		/**< untransposed matrix */
							const real_t* const M2,
							int_t dim0,					/**< leading dimension of untransposed M1 = leading dimension of M2 */
							int_t dim1,					/**< secondary dimension of untransposed M1 */
							int_t dim2,					/**< secondary dimension of M2 */
							boolean_t addToRes			/**< flag to specify whether to overwrite res, or simply add to it */
							);

/* ----------------------------------------------
 * Dense generic matrix-transposed matrix product
 * res = M1*M2.T
 *
 >>>>>                                            */
void multiplyMatrixMatrixTDenseDense( real_t* const res,
							const real_t* const M1,		/**< untransposed matrix */
							const real_t* const M2,
							int_t dim0,					/**< leading dimension of M1 */
							int_t dim1,					/**< secondary dimension of M1 = secondary dimension of untransposed M2 */
							int_t dim2					/**< leading dimension of untransposed M2 */
							);


/** Low-level scalar product */
real_t scalarProd(	const vector_t* const x,
					const vector_t* const y,
					int_t len 	);


return_t addVectors(	vector_t* const res,
						const vector_t* const x,
						const vector_t* const y,
						int_t len
						);


return_t addToVector(	vector_t* const res,
					const vector_t* const update,
					int_t len
					);


return_t subtractVectors(	vector_t* const res,
							const vector_t* const x,
							const vector_t* const y,
							int_t len
							);


return_t subtractFromVector(	vector_t* const res,
							const vector_t* const update,
							int_t len
							);


return_t negateVector(	vector_t* const res,
						int_t len
						);


return_t addMatrix(	matrix_t* const res,
					const matrix_t* const update,
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
real_t vectorNorm(	const vector_t* const vec,
					int_t len
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
return_t factorizePosDefMatrix( 	qpData_t* const qpData,
									matrix_t* const cholM,
									const matrix_t* const M,
									int_t dim0
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
return_t denseCholeskyFactorization(	qpData_t* const qpData,
										matrix_t* const cholM,
										const matrix_t* const M,
										int_t dim0
										);



return_t backsolveRT_ZTET(	qpData_t* const qpData,
							zx_matrix_t* const res,
							const zz_matrix_t* const RT,
							const zz_matrix_t* const ZT,
							x_vector_t* const sums,
							int_t dim0,						/* number of physical rows and columns in  RT (storage) = number of columns in ZT */
							int_t dim1						/* number of defined rows in ZT = number of defined rows and columns in RT */
							);

return_t backsolveRT_ZTCT(	qpData_t* const qpData,
							zx_matrix_t* const res,
							const zz_matrix_t* const RT,
							const zz_matrix_t* const ZTCT,
							x_vector_t* const sums,
							int_t dim0,						/* number of physical rows and columns in RT (storage) */
							int_t dim1						/* number of (well-defined) rows in ZTCT (same as ZT) */
							);


return_t addCInvHCT(	qpData_t* const qpData,
						xx_matrix_t* const res,
						const vv_matrix_t* const cholH,
						const xz_matrix_t* const C,
/*						const xx_matrix_t* const A,*/
/*						const xu_matrix_t* const B,*/
						const d2_vector_t* const y,
/*						xz_matrix_t* const C,*/			/**< temporary matrix to build up C as once */
						xx_matrix_t* const xxMatTmp,
						ux_matrix_t* const uxMatTmp,
						zx_matrix_t* const zxMatTmp
						);


return_t addMultiplyMatrixInvMatrixMatrixT(	qpData_t* const qpData,
											matrix_t* const res,
											const matrix_t* const cholM1,
											const matrix_t* const M2,
											const real_t* const y,			/**< vector containing non-zeros for columns of M2 to be eliminated */
											matrix_t* const Ztmp,			/**< temporary matrix of shape dim1 x dim0 */
											vector_t* const vecTmp,
											int_t dim0,						/**< dimensions of M2 */
											int_t dim1
											);



#endif	/* QP42_MATRIX_VECTOR_H */


/*
 *	end of file
 */
