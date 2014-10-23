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
 *	\file include/qp/dual_qp.h
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */


#ifndef QP42_DUAL_QP_H
#define QP42_DUAL_QP_H


#include <assert.h>
#include <qp/stage_qp_solver_clipping.h>
#include <qp/stage_qp_solver_qpoases.hpp>
#include <qp/types.h>
#include <qp/matrix_vector.h>
#include <qp/setup_qp.h>
#include <qp/qpdunes_utils.h>


return_t qpDUNES_solve(	qpData_t* const qpData
						);


void qpDUNES_logIteration( qpData_t* qpData,
						itLog_t* itLogPtr,
						real_t objValIncumbent,
						int_t lastActSetChangeIdx
						);



return_t qpDUNES_solveAllLocalQPs(	qpData_t* const qpData,
								const xn_vector_t* const lambda
								);


return_t qpDUNES_updateAllLocalQPs(	qpData_t* const qpData,
									const xn_vector_t* const lambda
									);


/* ----------------------------------------------
 * solve local QP
 * 
 >>>>>>                                           */
return_t qpDUNES_solveLocalQP(	qpData_t* const qpData,
							interval_t* const interval
							);


return_t qpDUNES_setupNewtonSystem(	qpData_t* const qpData
									);

return_t qpDUNES_factorNewtonSystem(	qpData_t* const qpData,
									boolean_t* const isHessianRegularized,
									int_t lastActSetChangeIdx
									);


return_t qpDUNES_computeNewtonGradient(	qpData_t* const qpData,
								xn_vector_t* gradient,
								x_vector_t* gradPiece
								);


return_t qpDUNES_factorizeNewtonHessian(	qpData_t* const qpData,
										xn2x_matrix_t* const cholHessian,
										xn2x_matrix_t* const hessian,
										boolean_t* isHessianRegularized
										);


return_t qpDUNES_factorizeNewtonHessianBottomUp(	qpData_t* const qpData,
												xn2x_matrix_t* const cholHessian,
												xn2x_matrix_t* const hessian,
												int_t lastActSetChangeIdx,
												boolean_t* isHessianRegularized
												);


return_t qpDUNES_solveNewtonEquation(	qpData_t* const qpData,
									xn_vector_t* const res,
									const xn2x_matrix_t* const cholHessian,	/**< lower triangular Newton Hessian factor */
									const xn_vector_t* const gradient
									);

return_t qpDUNES_solveNewtonEquationBottomUp(	qpData_t* const qpData,
											xn_vector_t* const res,
											const xn2x_matrix_t* const cholHessian,	/**< lower triangular Newton Hessian factor */
											const xn_vector_t* const gradient
											);

return_t qpDUNES_multiplyNewtonHessianVector(	qpData_t* const qpData,
											xn_vector_t* const res,
											const xn2x_matrix_t* const hessian, /**< Newton Hessian */
											const xn_vector_t* const vec	);


return_t qpDUNES_diffWorkingSet(	qpData_t* const qpData
								);


return_t qpDUNES_determineStepLength(	qpData_t* const qpData,
									xn_vector_t* const lambda,
									xn_vector_t* const deltaLambdaFS,
									uint_t* const itCntr,
									real_t* const alpha,
									real_t* const objValIncumbent,
									boolean_t newtonHessianRegularized
									);

return_t qpDUNES_backTrackingLineSearch(	qpData_t* const qpData,
										real_t* const alpha,
										uint_t* const itCntr,
/*										xn_vector_t* const lambda,*/
										const xn_vector_t* const deltaLambdaFS,
										xn_vector_t* const lambdaTry,
										int_t nV,
										real_t alphaMin,
										real_t alphaMax,
										real_t const objValIncumbent
										);

return_t qpDUNES_reductionLineSearchWithASChange(	qpData_t* const qpData,
													real_t* const alpha,
													uint_t* const itCntr,
													xn_vector_t* const lambda,
													const xn_vector_t* const deltaLambdaFS,
													xn_vector_t* const lambdaTry,
													int_t nV,
													real_t alphaMin,
													real_t alphaMax,
													real_t const objValIncumbent
													);

return_t qpDUNES_goldenSectionIntervalSearch(	qpData_t* const qpData,
											real_t* const alpha,
											uint_t* const itCntr,
											xn_vector_t* const lambda,
											const xn_vector_t* const deltaLambdaFS,
											xn_vector_t* const lambdaTry,
											int_t nV,
											real_t alphaMin,
											real_t alphaMax
											);

return_t qpDUNES_bisectionIntervalSearch(	qpData_t* const qpData,
										real_t* const alpha,
										uint_t* const itCntr,
/*										xn_vector_t* const lambda,*/
										const xn_vector_t* const deltaLambdaFS,
										xn_vector_t* const lambdaTry,
										int_t nV,
										real_t alphaMin,
										real_t alphaMax
										);

return_t qpDUNES_qpadApproxIntervalSearch(	qpData_t* const qpData,
										real_t* const alpha,
										uint_t* const itCntr,
										const xn_vector_t* const deltaLambdaFS,
										xn_vector_t* const lambdaTry,
										int_t nV,
										real_t alphaL,
										real_t alphaR
										);

return_t qpDUNES_gridSearch(	qpData_t* const qpData,
							real_t* const alpha,
							uint_t* const itCntr,
							real_t* const objValIncumbent,
							real_t alphaMin,
							real_t alphaMax
							);


return_t qpDUNES_infeasibilityCheck(	qpData_t* qpData
										);


void qpDUNES_getPrimalSol(	const qpData_t* const qpData,
							real_t* const z
							);

return_t qpDUNES_getDualSol(	const qpData_t* const qpData,
							real_t* const lambda,
							real_t* const y
							);


real_t qpDUNES_computeObjectiveValue(	qpData_t* const qpData
									);


real_t qpDUNES_computeParametricObjectiveValue(	qpData_t* const qpData,
												const real_t alpha
												);


uint_t qpDUNES_getActSet(	const qpData_t* const qpData,
						int_t *const *const actSetStatus
						);

uint_t qpDUNES_compareActSets(	const qpData_t* const qpData,
							const int_t *const *const newActSetStatus,
							const int_t *const *const oldActSetStatus,
							int_t *const lastActSetChangeIdx
							);


void qpDUNES_printIterationHeader( qpData_t* qpData );


void qpDUNES_printIteration( 	qpData_t* qpData,
							itLog_t* itLogPtr
							);


#endif	/* QP42_DUAL_QP_H */


/*
 *	end of file
 */
