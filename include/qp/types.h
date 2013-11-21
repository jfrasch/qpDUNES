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
 *	\file include/qp/types.h
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 *
 *	Declaration of all non-built-in types.
 */


#ifndef QPDUNES_TYPES_H
#define QPDUNES_TYPES_H


#ifdef __MATLAB__
#include "matrix.h"	/* for mwSize types */
#endif


/*#define __USE_ASSERTS__*/					/* check return values of routines rigorously (useful for debugging) */
/*#undef __USE_ASSERTS__*/

/*#define __SUPPRESS_ALL_OUTPUT__*/			/* no printing */
/*#undef __SUPPRESS_ALL_OUTPUT__*/

/*#define __SUPPRESS_ALL_WARNINGS__*/		/* do not display warnings */
/*#undef __SUPPRESS_ALL_WARNINGS__*/

/*#define __MEASURE_TIMINGS__*/				/* measure computation times */
/*#undef __MEASURE_TIMINGS__*/

#define __ANALYZE_FACTORIZATION__			/* log inverse Newton Hessian for analysis */
#undef __ANALYZE_FACTORIZATION__

#define __QPDUNES_PARALLEL__				/* use openMP parallelization */
#undef __QPDUNES_PARALLEL__

#define PRINTING_PRECISION 14

#ifdef __MATLAB__
	#define MAX_STR_LEN 2560
#endif


/** colored/formatted terminal output */
#ifndef __MATLAB__
	#define COL_STD  "\033[0m"
	#define COL_SUCC "\033[92m"
	#define COL_WARN "\033[93m"
	#define COL_ERR  "\033[91m"
#else
	#define COL_STD  " ***"
	#define COL_SUCC "*** "
	#define COL_WARN "*** "
	#define COL_ERR  "*** "
#endif


/** simple types */
#ifndef __MATLAB__
	#ifndef int_t
		typedef int int_t;
	#endif
	#ifndef uint_t
		typedef unsigned int uint_t;
	#endif
	#ifndef real_t
		#ifdef __USE_SINGLE_PRECISION__
			typedef float real_t;
		#else
			typedef double real_t;
		#endif	/* __USE_SINGLE_PRECISION__ */
	#endif
#else	/* __MATLAB__ */
	typedef int int_t;
	typedef unsigned int uint_t;
	typedef double real_t;
#endif	/* __MATLAB__ */


#if !defined(__STATIC_MEMORY__)
	#define _NX_ (qpData->nX)
	#define _NU_ (qpData->nU)
	#define _NZ_ (qpData->nZ)
	#define _NV( I ) (qpData->intervals[ I ]->nV)
	#define _NI_ (qpData->nI)
	#define _ND( I ) (qpData->intervals[ I ]->nD)
	#define _NDTTL_ (qpData->nDttl)
#endif



/** MATRIX ACCESS */
/*                                                block offset   row offset   column offset (0=diag,-1=supDiag)   column */
#define accHessian( K, L, I, J )	hessian->data[ (K)*2*_NX_*_NX_  + (I)*2*_NX_   + (1+L)*_NX_                         + (J) ]
/*                                                        block offset   row offset   column offset (0=diag,-1=supDiag)   column */
#define accCholHessian( K, L, I, J )	cholHessian->data[ (K)*2*_NX_*_NX_  + (I)*2*_NX_   + (1+L)*_NX_                         + (J) ]

#define accH( I, J )	H->data[ (I)*nV + (J) ]

#define accC( I, J )	C->data[ (I)*_NZ_ + (J) ]

#define accM( I, J, DIM )	M[ (I)*(DIM) + (J) ]		/**< generic low level matrix access */
#define accMT( I, J, DIM )	M[ (J)*(DIM) + (I) ]		/**< generic low level transposed matrix access */
#define accL( I, J, DIM )	L[ (I)*(DIM) + (J) ]		/**< generic low level lower triangular matrix access */


/** Advanced types */

/** A simple Boolean type */
#ifdef __APPLE__
	typedef unsigned int boolean_t;
	#ifndef QPDUNES_FALSE
		#define QPDUNES_FALSE 0
	#endif
	#ifndef QPDUNES_TRUE
		#define QPDUNES_TRUE 1
	#endif
#else
	typedef enum
	{
		QPDUNES_FALSE,					/**< ... */
		QPDUNES_TRUE					/**< ... */
	} boolean_t;
#endif


/** Matrix sparsity type */
typedef enum
{
	QPDUNES_MATRIX_UNDEFINED,		/**< ... */
	QPDUNES_DENSE,					/**< ... */
	QPDUNES_SPARSE,					/**< ... */
	QPDUNES_DIAGONAL,				/**< ... */
	QPDUNES_IDENTITY,				/**< ... */
	QPDUNES_ALLZEROS					/**< ... */
} sparsityType_t;


/** Log level */
typedef enum
{
	QPDUNES_LOG_OFF = 0,					/**< ... */
	QPDUNES_LOG_ITERATIONS,					/**< ... */
	QPDUNES_LOG_ALL_DATA					/**< ... */
} logLevel_t;


/** Local QP solver list */
typedef enum
{
	QPDUNES_STAGE_QP_SOLVER_UNDEFINED,		/**< ... */
	QPDUNES_STAGE_QP_SOLVER_CLIPPING,		/**< ... */
	QPDUNES_STAGE_QP_SOLVER_QPOASES			/**< ... */
} qp_solver_t;


/** Newton Hessian regularization types */
typedef enum
{
	QPDUNES_REG_LEVENBERG_MARQUARDT,				/**< 0 = ... */
	QPDUNES_REG_NORMALIZED_LEVENBERG_MARQUARDT,		/**< 1 = ... */
	QPDUNES_REG_SINGULAR_DIRECTIONS,				/**< 2 = regularize only in singular directions during Cholesky factorization */
	QPDUNES_REG_UNCONSTRAINED_HESSIAN,				/**< 3 = ... */
	QPDUNES_REG_GRADIENT_STEP						/**< 4 = ... */
} nwtnHssnRegType_t;


/** Newton Hessian factorization  algorithms */
typedef enum
{
	QPDUNES_NH_FAC_BAND_FORWARD,		/**< 0 = ... */
	QPDUNES_NH_FAC_BAND_REVERSE			/**< 1 = ... */
} nwtnHssnFacAlg_t;


/** Line search types */
typedef enum
{
	QPDUNES_LS_BACKTRACKING_LS,						/**< 0 = ... */
	QPDUNES_LS_BACKTRACKING_LS_WITH_AS_CHANGE,		/**< 1 = fast backtracking until progress, bisection until AS change */
	QPDUNES_LS_GOLDEN_SECTION_LS,					/**< 2 = ... */
	QPDUNES_LS_GRADIENT_BISECTION_LS,				/**< 3 = ... */
	QPDUNES_LS_ACCELERATED_GRADIENT_BISECTION_LS,	/**< 4 = fast backtracking first, then gradient based bisection for refinement */
	QPDUNES_LS_GRID_LS,								/**< 5 = evaluate objective function on a grid and take minimum */
	QPDUNES_LS_ACCELERATED_GRID_LS					/**< 6 = fast backtracking first, then grid search for refinement */
} lineSearchType_t;


/** Error codes */
typedef enum
{
	QPDUNES_UNTERMINATED,
	QPDUNES_OK = 0,							/**< ... */

	QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND,
	QPDUNES_ERR_STAGE_QP_INFEASIBLE,
	QPDUNES_ERROR_STAGE_COUPLING_INFEASIBLE,

	QPDUNES_ERR_UNKNOWN_ERROR,
	QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE,	/**< ... */
	QPDUNES_ERR_UNKNOWN_LS_TYPE,
	QPDUNES_ERR_INVALID_ARGUMENT,
	QPDUNES_ERR_ITERATION_LIMIT_REACHED,
	QPDUNES_ERR_DIVISION_BY_ZERO,
	QPDUNES_ERR_NUMBER_OF_MAX_LINESEARCH_ITERATIONS_REACHED,
	QPDUNES_ERR_DECEEDED_MIN_LINESEARCH_STEPSIZE,
	QPDUNES_ERR_EXCEEDED_MAX_LINESEARCH_STEPSIZE,
	QPDUNES_ERR_NEWTON_SYSTEM_NO_ASCENT_DIRECTION,
	QPDUNES_NOTICE_NEWTON_MATRIX_NOT_SET_UP
} return_t;



/** 
 *	\brief generic matrix data type
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
typedef struct
{
	/** matrix property flags */
	sparsityType_t sparsityType;
/* 	boolean_t isDefined;
 	boolean_t hasChanged;*/

	/** matrix data array */
	real_t* data;
} matrix_t;

typedef matrix_t xx_matrix_t;
typedef matrix_t xu_matrix_t;
typedef matrix_t xz_matrix_t;
typedef matrix_t ux_matrix_t;
typedef matrix_t uu_matrix_t;
typedef matrix_t zx_matrix_t;
typedef matrix_t zz_matrix_t;
typedef matrix_t vv_matrix_t;
typedef matrix_t dz_matrix_t;

/**
 * Special Newton hessian storage format:
 *
 *  A symmetric block tri-diagonal Newton Hessian
 * 		[ D L'       ]
 * 		[ L D L'     ]
 * 		[   L D L'   ]
 * 		[      ...   ]
 * 		[        L D ]
 *
 *  is stored as
 * 		[ - D ]
 * 		[ L D ]
 * 		[ L D ]
 * 		[ ... ]
 * 		[ L D ]
 *
 */
typedef matrix_t xn2x_matrix_t;

typedef matrix_t xnxn_matrix_t;



/** 
 *	\brief generic vector data type
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
typedef struct
{
	/** vector property flags */
	boolean_t isDefined;
	boolean_t hasChanged;

	/** vector data array */
	real_t* data;
} vector_t;

typedef vector_t x_vector_t;
typedef vector_t u_vector_t;
typedef vector_t z_vector_t;
typedef vector_t d_vector_t;
typedef vector_t d2_vector_t;
typedef vector_t d2n1_vector_t;
typedef vector_t xn_vector_t;
typedef vector_t zn_vector_t;
typedef vector_t zn1_vector_t;



/** 
 *	\brief generic integer vector data type
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
typedef struct
{
	/** vector property flags */
	boolean_t isDefined;
	boolean_t hasChanged;

	/** vector data array */
	int_t* data;
} intVector_t;

typedef intVector_t zn_intVector_t;



/**
 *	\brief pointer to qpOASES object for C++ method access
 *
 *	...
 *
 *	\author Janick Frasch
 *	\version 1.0beta
 *	\date 2013
 */
typedef void qpoasesObject_t;


/**
 *	\brief struct with auxiliary data for QPOASES QP solver
 *
 *	...
 *
 *	\author Janick Frasch
 *	\version 1.0beta
 *	\date 2013
 */
typedef struct
{
	qpoasesObject_t* qpoasesObject;

	/* workspace */
	z_vector_t qFullStep;			/**< linear term corresponding to full-step in lambda */
	real_t pFullStep;				/**< constant term corresponding to full-step in lambda */
} qpSolverQpoases_t;



/**
 *	\brief struct with auxiliary data for clipping QP solver
 *
 *	...
 *
 *	\author Janick Frasch
 *	\version 1.0beta
 *	\date 2013
 */
typedef struct
{
	z_vector_t zUnconstrained;	/**< unconstrained primal solution for current lambda guess */
	z_vector_t dz;				/**< delta z - update in primal variables corresponding to a full step deltaLambda */

	/* workspace */
	z_vector_t qStep;			/**< step in linear term for line search */
	real_t pStep;				/**< step in constant term for line search */
} qpSolverClipping_t;


/**
 *	\brief Hessian interval data type and dynamic constraint interval data type
 *
 *	Datatype for one interval of the QP data.
 * 
 *  Structure of the Hessian block:
 *   (x)' (Q  S) (x)
 *   (u)  (S' R) (u)
 *
 *  Structure of interval dynamics:
 *   x_{k+1} = A_k*x_k + B_k*u_k + c_k
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
typedef struct
{
	uint_t id;					/**< stage index */


	/* dimensions */
	uint_t nD;					/**< number of constraints */
	uint_t nV;					/**< number of variables */
	

	/* primal objective function */
	vv_matrix_t H;				/**< Hessian */
	vv_matrix_t cholH;			/**< inverse of Hessian */
	real_t HQNorm;				/**< norm of Q-part of Hessian */     /*FIXME: choose which norm to compute exactly, etc.*/
	
	z_vector_t g;				/**< primal gradient block */
	

	/* dualized objective function */
	z_vector_t q;				/**< linear objective function term after dualization */
	real_t p;					/**< constant objective function term after dualization */


	/* dynamic system */
	xz_matrix_t C;				/**< u-part of constraint matrix */
	x_vector_t  c;				/**< constant part */


	/* constraints */
	z_vector_t  zLow;			/**< lower variable bound */
	z_vector_t  zUpp;			/**< upper variable bound */
	dz_matrix_t D;				/**< full constraint matrix */
	d_vector_t  dLow;			/**< constraint lower bound */
	d_vector_t  dUpp;			/**< constraint upper bound */
	

	/* primal QP solution */
	z_vector_t z;				/**< full primal solution for current lambda guess */
	real_t optObjVal;			/**< objective value */


	/* dual QP solution */
	d2_vector_t y;				/**< stage constraint multiplier vector  */


	/* QP solver */
	qp_solver_t qpSolverSpecification;		/**< dedicated QP solver */
	
	qpSolverClipping_t qpSolverClipping;	/**< workspace for clipping QP solver */
	qpSolverQpoases_t qpSolverQpoases;		/**< pointer to qpOASES object */
	
	boolean_t actSetHasChanged;				/**< indicator flag whether an active set change occurred on this
										     	 interval during the current iteration */


	/* workspace */
	x_vector_t lambdaK;		/**<  */
	x_vector_t lambdaK1;	/**<  */
	
	x_vector_t xVecTmp;			/**<  */
	u_vector_t uVecTmp;			/**<  */
	z_vector_t zVecTmp;			/**<  */
	
} interval_t;



/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
typedef struct
{
	/* iteration limits */
	int_t maxIter;
	int_t maxNumLineSearchIterations;			/**< maximum number of line search steps in solution of Newton system */
	int_t maxNumLineSearchRefinementIterations;	/**< maximum number of refinement line search steps to find point with AS change */

	/* printing */
	int_t printLevel;							/**< Amount of information printed:   0 = no output
																					  1 = only errors and success
																					  2 = additionally iterations and warnings
																					  3 = debug information
																					  */
	/* logging */
	logLevel_t logLevel;						/**< Amount of information logged */

	int_t printIntervalHeader;
	boolean_t printIterationTiming;
	boolean_t printLineSearchTiming;
	/* TODO: use print precision here for variable precision matrix and vector printing */

	/* numerical tolerances */
	real_t stationarityTolerance;
	real_t equalityTolerance;
	real_t newtonHessDiagRegTolerance;	/**< Tolerance on diagonal elements of Newton Hessian before regularizing */
	real_t activenessTolerance;			/**< Tolerance in constraint violation before a constraint is considered active => needed to avoid weakly active constraints, which might yield suboptimal Hessian directions */
	
	real_t QPDUNES_ZERO;					/** Numerical value of zero (for situations in which it would be unreasonable to compare with 0.0). Has to be positive. */
	real_t QPDUNES_INFTY;					/**< Numerical value of infinity (e.g. for non-existing bounds). Has to be positive. */
	real_t ascentCurvatureTolerance;	/**< Tolerance when a step is called a zero curvature step */
	
	/* additional options */
	int_t nbrInitialGradientSteps;		/**< after the first Newton step a number of cheaper gradient
											 steps with line search can be used to drive the method
											 faster to the solution */
	boolean_t checkForInfeasibility;	/**< perform checks for infeasibility of the problem */

	/* regularization options */
	nwtnHssnRegType_t regType;
	real_t regParam;					/**< Levenberg-Marquardt relaxation parameter */
	
	nwtnHssnFacAlg_t nwtnHssnFacAlg;

	/* line search options */
	lineSearchType_t lsType;
	real_t lineSearchReductionFactor;
	real_t lineSearchIncreaseFactor;
	real_t lineSearchMinAbsProgress;
	real_t lineSearchMinRelProgress;
	real_t lineSearchStationarityTolerance;
	real_t lineSearchMaxStepSize;
	int_t lineSearchNbrGridPoints;		/**< number of grid points for grid line search */

	/* qpOASES options */
	real_t qpOASES_terminationTolerance;

} qpOptions_t;



/**
 *	\brief log type for single iteration
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
typedef struct
{
	/* data: vectors, matrices */
	xn_vector_t lambda;
	xn_vector_t deltaLambda;

	xn2x_matrix_t hessian;
	xn2x_matrix_t cholHessian;
	xn_vector_t gradient;

	zn1_vector_t dz;
	zn1_vector_t zUnconstrained;
	zn1_vector_t z;

	d2n1_vector_t y;

	xn_vector_t regDirections; 		/* log the regularized directions delta, for Newton system (H+diag(delta))*lambda = -g */

	#ifdef __ANALYZE_FACTORIZATION__
	xnxn_matrix_t invHessian;
	#endif

	int_t** ieqStatus;
	int_t** prevIeqStatus;

	/* timings */
	real_t tIt;
	real_t tNwtnSetup;
	real_t tNwtnSolve;
	real_t tQP;
	real_t tLineSearch;

	/* statuses */
	real_t gradNorm;
	real_t lambdaNorm;
	real_t stepNorm;
	real_t stepSize;
	real_t objVal;

	uint_t nActConstr;
	uint_t nChgdConstr;
	int_t lastActSetChangeIdx;


	/* flags, etc. */
	uint_t itNbr;
	boolean_t isHessianRegularized;
	uint_t numLineSearchIter;

} itLog_t;


/**
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
typedef struct
{
	/* Problem dimensions */
	uint_t nI;
	uint_t nX;
	uint_t nU;
	uint_t nZ;
	uint_t nDttl;				/**< total number of local constraints */

	/* Problem data */
	interval_t** intervals;

	/* options */
	qpOptions_t qpOptions;

	/* iterations log */
	itLog_t* itLog;

	int_t numIter;

} log_t;



/** 
 *	\brief ...
 *
 *	...
 *
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 */
typedef struct
{
	/* variables */
	uint_t nI;
	uint_t nX;
	uint_t nU;
	uint_t nZ;
	uint_t nDttl;				/**< total number of local constraints */

	interval_t** intervals;		/**< array of pointers to interval structs; double pointer for more efficient shifting */

	xn_vector_t lambda;
	xn_vector_t deltaLambda;
	
	xn2x_matrix_t hessian;
	xn2x_matrix_t cholHessian;
	xn_vector_t gradient;
	

	real_t alpha;
	real_t optObjVal;
	
	qpOptions_t options;
	
	
	/* workspace */
	x_vector_t xVecTmp;			/**<  */
	u_vector_t uVecTmp;			/**<  */
	z_vector_t zVecTmp;			/**<  */
	xn_vector_t xnVecTmp;		/**<  */
	xn_vector_t xnVecTmp2;		/**<  */
	
	xx_matrix_t xxMatTmp;		/**<  */
	xx_matrix_t xxMatTmp2;		/**<  */
	ux_matrix_t uxMatTmp;		/**<  */
	xz_matrix_t xzMatTmp;		/**<  */
	zx_matrix_t zxMatTmp;		/**<  */
	zz_matrix_t zzMatTmp;		/**<  */
	zz_matrix_t zzMatTmp2;		/**<  */
	
	/* log */
	log_t log;

} qpData_t;


#endif	/* QPDUNES_TYPES_H */


/*
 *	end of file
 */
