/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.1 beta
|
|   ./src/SolverCG.cpp
|
|                     Written by T.Takeda,    2012/03/26
|                                Y.Sato,      2012/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include "SolverCG.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
#include <cmath>
#ifdef _OPENMP
#include "omp.h"
#endif

namespace pmw
{
CSolverCG::CSolverCG(iint iter_max = 100,
		double tolerance = 1.0e-8,
		iint method = 1,
		iint precondition = 1,
		bool flag_iter_log = false,
		bool flag_time_log = false)
	: CSolver(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log)
{
}
CSolverCG::~CSolverCG()
{
}
////uiint CSolverCG::doSolve_(const CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
////}
uiint CSolverCG::doSolve( CAssyMatrix *pA, const CAssyVector *pB, CAssyVector *pX,
		          iint iter_max, double tolerance, bool flag_iter_log, bool flag_time_log)
{
	CHecMPI *pMPI= CHecMPI::Instance();
	uiint myrank = pMPI->getRank();
	/*
	if(myrank == 0){
		pA->dump();
		pB->dump();
	}
	*/
        Utility::CLogger *pLogger= Utility::CLogger::Instance();

	double bnrm2_inv, bnrm2, rnrm2;
	double alpha, beta, rho, rho_prev, resid;
	iint iter, iter_precond; 
#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CSolverCG::doSolve1111 \n");
#endif
   	printf(" --- start of CG solver --- %ld \n",getPrecondition());
	uiint len = pB->size();
	CAssyVector R(pX); 
	CAssyVector Z(pX); 
	CAssyVector P(pX); 
	CAssyVector &Q = Z;

	pA->residual(pX, pB, &R);//////////////////////////// 内部:AssyMatrix::multVector

        bnrm2= pB->norm2();

        if(bnrm2==0.0){
            bnrm2_inv= 0.0;      //----( 右辺ベクトル== 0.0 )
        }else{
            bnrm2_inv= 1.0/bnrm2;//----( 右辺ベクトル > 0.0 )
        }
        
	iint itype = getPrecondition();// <-- mPrecondition
	pA->setupPreconditioner( itype );

	Z.setZero();
        
        setCGridSolver( pA, iter_max, tolerance, flag_iter_log, flag_time_log);//コースグリッド線形ソルバー生成・セット

	for (iter = 0; iter < getIterMax(); iter++) {
            /*
            if(myrank == 0){
                    pX->dump();
            }
            */

            iter_precond = 1;
            if( itype == 2 ) {
                Z.subst(pX);
                //---- R:右辺、Z:解
                pA->MGInitialGuess(&R, &Z);//------ MultiGrid前処理
                //----
            } else {
                pA->precond(&R, &Z, iter_precond);
            }

            // "nan"判定,  Z:残差の解
            #ifdef MSVC
            if(_isnan(Z.norm2()) )
            #else
            if( isnan(Z.norm2()) )
            #endif
                pLogger->Info_format(Utility::LoggerMode::Error, "%s", "SolverCG,  Z.norm2 is NaN");

            rho = R.innerProd(&Z);
            if (iter > 0) {
                beta = rho / rho_prev;
                #pragma omp parallel for
                for(uiint i = 0; i < len; i++)
                    P[i] = Z[i] + beta * P[i];
            } else {
                P.subst(&Z);
            }
            /////////////////////////////////////////////////////////////
            pA->multVector(&P, &Q);//--------- 行列・ベクトル積(通信&MPC接合)
            /////////////////////////////////////////////////////////////

            alpha = rho / P.innerProd(&Q);
            
            // "nan"判定: alpha
            #ifdef MSVC
            if(_isnan(alpha) ){
            #else
            if( isnan(alpha) ){
            #endif
                pLogger->Info_format(Utility::LoggerMode::Error, "%s", "SolverCG,  alpha is NaN");
                alpha= 1.0;
            }
            
            #pragma omp parallel for
            for (uiint i = 0; i < len; i++) {
                (*pX)[i] += alpha * P[i];//--------------------- X = X + α*P
                R[i]     -= alpha * Q[i];//--------------------- 残差
            }

            rnrm2 = R.norm2();//-------- 残差norm2
            if(rnrm2==0){
                resid= 0.0;
            }else{
                resid = sqrt(rnrm2 * bnrm2_inv);//---残差
            }
            
            if(myrank == 0){
            	if (getFlagIterLog()) {
                printf("%5ld %16.6e\n", iter + 1, resid);
            	}
            	printf("iteration:%5ld, residue:%16.6e \n", iter + 1, resid);
            }
            rho_prev = rho;

            if (resid < getTolerance()) break;
            
	};//iter end

////        //----- この処理は、本来ここ => MPCデバッグのため後ろのまま.
////    P.subst(pX);
////        ///////////////////////////////////////////////////////////
////	pA->multMPC(&P, pX);//------------------------ u = Tu' + ug
////        ///////////////////////////////////////////////////////////

	if( iter == getIterMax() ){
            cout << "---------------------- IterMax iter:" << iter << " " << endl;
            deleteCGridSolver();//コースグリッドソルバー削除
            return 0;
        }

	P.subst(pX);
        ///////////////////////////////////////////////////////////
	pA->multMPC(&P, pX);//------------------------ u = Tu' + ug
        ///////////////////////////////////////////////////////////

   	printf(" --- end of CG solver --- \n");
        
        deleteCGridSolver();//コースグリッドソルバー削除

#ifdef ADVANCESOFT_DEBUG
	printf(" ----- the solution vector of solver ----- \n");
	for(uiint i=0; i < 5; i++)
		printf(" %d : %e %e %e \n", i, (*pX)[i][0], (*pX)[i][1], (*pX)[i][2]);
	for(uiint i=len-5; i < len; i++)
		printf(" %d : %e %e %e \n", i, (*pX)[i][0], (*pX)[i][1], (*pX)[i][2]);
   	printf(" exit CSolverCG::doSolve \n");
#endif
	return 1;
}
}
