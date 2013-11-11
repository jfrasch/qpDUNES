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
 *	\file src/stage_qp_solver_qpoases.c
 *	\author Janick Frasch
 *	\version 1.0beta
 *	\date 2012
 */


#include <qp/stage_qp_solver_qpoases.hpp>


extern "C" {

/* ----------------------------------------------
 * solve QP
 *
#>>>>>>                                           */
qpoasesObject_t* qpOASES_contructor( qpData_t* qpData,
									 int_t nV,
							 	 	 int_t nC )
{

}
/*<<< END OF qpOASES_contructor */


/* ----------------------------------------------
 * solve QP
 *
#>>>>>>                                           */
void qpOASES_destructor( qpoasesObject_t** qpoasesObject )
{

}
/*<<< END OF qpOASES_destructor */


/* ----------------------------------------------
 * first QP solution
 *
#>>>>>>                                           */
return_t qpOASES_setup( qpData_t* qpData,
						qpoasesObject_t* qpoasesObject,
						interval_t* interval,
						zz_matrix_t* H,
						z_vector_t* g,
						d_vector_t* zLow,
						d_vector_t* zUpp,
						dz_matrix_t* D,
						d_vector_t* dLow,
						d_vector_t* dUpp	)
{

}
/*<<< END OF qpOASES_setup */



/* ----------------------------------------------
 * update QP data
 *
 * qStep = C.T*lambdaK1 - [lambdaK.T 0]
 * pStep = c*lambdaK1
#>>>>>>                                           */
return_t qpOASES_updateStageData(	qpData_t* const qpData,
									interval_t* const interval,
									const z_vector_t* const lambdaK,
									const z_vector_t* const lambdaK1
									)
{

}
/*<<< END OF qpOASES_updateStageData */


/* ----------------------------------------------
 * Hotstarted QP solution for updated gradient
 *
#>>>>>>                                           */
return_t qpOASES_hotstart( 	qpData_t* qpData,
							qpoasesObject_t* qpoasesObject,
							interval_t* interval,
							z_vector_t* q 				)
{

}
/*<<< END OF qpOASES_hotstart */


/* ----------------------------------------------
 * Get projected Hessian
 *
#>>>>>>                                           */
return_t qpOASES_getCholZTHZ( 	qpData_t* qpData,
								qpoasesObject_t* qpoasesObject,
								zz_matrix_t* cholZTHZ 				)
{

}
/*<<< END OF qp42_solveLocalQP */


/* ----------------------------------------------
 * Get null-space basis matrix Z
 *
#>>>>>>                                           */
return_t qpOASES_getZT( qpData_t* qpData,
						qpoasesObject_t* qpoasesObject,
						int_t* nFree,
						zz_matrix_t* ZT 				)
{

}
/*<<< END OF qp42_solveLocalQP */


/* ----------------------------------------------
 * do a step of length alpha
 *
#>>>>>>                                           */
return_t qpOASES_doStep( qpData_t* const qpData,
						 qpoasesObject_t* qpoasesObject,
						 interval_t* const interval,
						 real_t alpha,
						 z_vector_t* const z,
						 d2_vector_t* const mu,
						 z_vector_t* const qCandidate,
						 real_t* const p				)
{

}
/*<<< END OF qpOASES_doStep */




} /* extern C */

/*
 *	end of file
 */
