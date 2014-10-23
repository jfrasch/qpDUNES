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
 *	\file include/qp/stage_qp_solver_qpoases.h
 *	\author Janick Frasch
 *	\version 1.0beta
 *	\date 2012
 */


#ifndef QP42_STAGE_QP_SOLVER_QPOASES_H
#define QP42_STAGE_QP_SOLVER_QPOASES_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */


#include <qp/types.h>
#include <qp/matrix_vector.h>
#include <qp/qpdunes_utils.h>


/** ... */
qpoasesObject_t* qpOASES_contructor(	qpData_t* qpData,
										int_t nV,
										int_t nC );


/** ... */
void qpOASES_destructor( qpoasesObject_t** qpoasesObject );


/** ... */
return_t qpOASES_setup( qpData_t* qpData,
						qpoasesObject_t* qpoasesObject,
						interval_t* interval,
						zz_matrix_t* H,
						z_vector_t* g,
						d_vector_t* zLow,
						d_vector_t* zUpp,
						dz_matrix_t* D,
						d_vector_t* dLow,
						d_vector_t* dUpp	);


/** ... */
return_t qpOASES_updateStageData(	qpData_t* const qpData,
									interval_t* const interval,
									const z_vector_t* const lambdaK,
									const z_vector_t* const lambdaK1
									);


/** ... */
return_t qpOASES_hotstart( 	qpData_t* qpData,
							qpoasesObject_t* qpoasesObject,
							interval_t* interval,
							z_vector_t* q 				);


/** ... */
return_t qpOASES_getCholZTHZ( 	qpData_t* qpData,
								qpoasesObject_t* qpoasesObject,
								zz_matrix_t* cholZTHZ 				);


/* ----------------------------------------------
 * Get null-space basis matrix Z
 *
#>>>>>>                                           */
return_t qpOASES_getZT( qpData_t* qpData,
						qpoasesObject_t* qpoasesObject,
						int_t* nFree,
						zz_matrix_t* ZT 				);


return_t qpOASES_getPrimalDualVariables( qpData_t* qpData,
										 qpoasesObject_t* qpoasesObject 	);


/** ... */
return_t qpOASES_doStep( qpData_t* const qpData,
						 qpoasesObject_t* qpoasesObject,
						 interval_t* const interval,
						 real_t alpha,
						 z_vector_t* const z,
						 d2_vector_t* const mu,
						 z_vector_t* const q,
						 real_t* const p				);


#ifdef __cplusplus
} /* extern C */
#endif /* __cplusplus */

#endif	/* QP42_STAGE_QP_SOLVER_QPOASES_H */


/*
 *	end of file
 */
