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
 *	\file include/qp/stage_qp_solver_clipping.h
 *	\author Janick Frasch
 *	\version 1.0beta
 *	\date 2012
 */


#ifndef QP42_STAGE_QP_SOLVER_CLIPPING_H
#define QP42_STAGE_QP_SOLVER_CLIPPING_H


#include <qp/types.h>
#include <qp/matrix_vector.h>
#include <qp/qpdunes_utils.h>
#include <assert.h>

/** ... */
return_t directQpSolver_solve( qpData_t* qpData, 
									 interval_t* interval );

/** ... */
return_t clippingQpSolver_updateStageData(	qpData_t* const qpData,
												interval_t* const interval,
												const z_vector_t* const lambdaK,
												const z_vector_t* const lambdaK1
												);


/** ... */
return_t directQpSolver_solveUnconstrained( qpData_t* const qpData,
									 	 	 	 interval_t* const interval,
									 	 	 	 const z_vector_t* const qStep );


/** ... */
return_t directQpSolver_getMinStepsize( const qpData_t* const qpData,
										const interval_t* const interval,
										real_t* alphaMin );


/** ... */
return_t directQpSolver_doStep( qpData_t* const qpData,
								interval_t* const interval,
								const z_vector_t* const stepDir,
								real_t alpha,
								z_vector_t* const zUnconstrained,
								z_vector_t* const z,
								d2_vector_t* const mu,
								z_vector_t* const q,
								real_t* const p				);



/** ... */
return_t directQpSolver_saturateVector(	qpData_t* const qpData,
										d_vector_t* const vec, 
										d2_vector_t* const mu,
										const d_vector_t* const lb,
										const d_vector_t* const ub,
										const zz_matrix_t* const H,
										int_t nD
										);


/** ... */
return_t clippingQpSolver_ratioTest(	qpData_t* const qpData,
									real_t* minStepSizeASChange,	/* minimum step size that leads to active set change */
									d_vector_t* const zStepDir,
									d2_vector_t* const mu,			/* pseudo multipliers, resembling the gaps to the bounds; + active, - inactive */
									const d_vector_t* const lb,
									const d_vector_t* const ub,
									int_t nV
									);


/** ... */
real_t directQpSolver_getObjectiveValue(	qpData_t* const qpData,
												interval_t* const interval
												);


#endif	/* QP42_STAGE_QP_SOLVER_CLIPPING_H */


/*
 *	end of file
 */
