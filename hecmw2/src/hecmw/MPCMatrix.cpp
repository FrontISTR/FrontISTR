/*
 * MPCMatrix.cpp
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
 */

#include "MPCMatrix.h"
#include "AssyVector.h"
#include "Equation.h"

namespace pmw
{

CMPCMatrix::CMPCMatrix() :
	mnEquation(0)
{
	// TODO Auto-generated constructor stub
}

CMPCMatrix::~CMPCMatrix()
{
	// TODO Auto-generated destructor stub
}

typedef std::vector<CEquation*> CVEqn;
typedef CVEqn::iterator CVEqnIter;
typedef CVEqn::const_iterator CVEqnConstIter;

typedef std::vector<CEqnTerm> CVEqTrm;
typedef CVEqTrm::iterator CVEqTrmIter;
typedef CVEqTrm::const_iterator CVEqTrmConstIter;

void CMPCMatrix::multVector(CAssyVector *pV, CAssyVector *pP) const
{
    pP->subst(pV);

    for (CVEqnConstIter i_eqn = mvEquation.begin(); i_eqn != mvEquation.end(); i_eqn++) {
        const CVEqTrm &vTerm = (*i_eqn)->getTerm();

        CVEqTrmConstIter j_term = vTerm.begin();
        double &slave_term = (*pP)(j_term->meshID(), j_term->nodeID(), j_term->dof());
        j_term++;

        double sum= 0.0;
        for (; j_term != vTerm.end(); j_term++) {
                sum += j_term->coef() * (*pV)(j_term->meshID(), j_term->nodeID(), j_term->dof());
        }
        slave_term = sum;
    }
}

void CMPCMatrix::transMultVector(CAssyVector *pV, CAssyVector *pP) const
{
    pP->subst(pV);

    for(CVEqnConstIter i_eqn = mvEquation.begin(); i_eqn != mvEquation.end(); i_eqn++){
        const CVEqTrm &vTerm = (*i_eqn)->getTerm();

        CVEqTrmConstIter j_term = vTerm.begin();
        (*pP)(j_term->meshID(), j_term->nodeID(), j_term->dof()) = 0.0;
        double slave_val = (*pV)(j_term->meshID(), j_term->nodeID(), j_term->dof());
        j_term++;

        for(; j_term != vTerm.end(); j_term++){
            (*pP)(j_term->meshID(), j_term->nodeID(), j_term->dof()) += j_term->coef() * slave_val;
        };
    };
}

// デバッグ
// ------
// 2011.01.13
//
void CMPCMatrix::dump()
{
    uiint nNumOfEqu = mvEquation.size();
    uiint i,ii;
    for(i=0; i < nNumOfEqu; i++){
        uiint nNumOfTerm;
        nNumOfTerm = mvEquation[i]->numTerm();
        
        cout << " " ;
        for(ii=0; ii < nNumOfTerm; ii++){
            cout << " NodeID=" << mvEquation[i]->getTerm(ii).nodeID();
            cout << " DOF=" << mvEquation[i]->getTerm(ii).dof();
            cout << " Coef=" << mvEquation[i]->getTerm(ii).coef();
        };
        cout << endl;
    };
}

}//namespace pmw
