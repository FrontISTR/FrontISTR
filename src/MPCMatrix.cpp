/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   MPCMatrix.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
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
}
CMPCMatrix::~CMPCMatrix()
{
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
		double sum = 0.0;
		for (; j_term != vTerm.end(); j_term++) {
			sum += j_term->coef() * (*pV)(j_term->meshID(), j_term->nodeID(), j_term->dof());
		}
		slave_term = sum;
	}
}
void CMPCMatrix::transMultVector(CAssyVector *pV, CAssyVector *pP) const
{
	pP->subst(pV);
	for (CVEqnConstIter i_eqn = mvEquation.begin(); i_eqn != mvEquation.end(); i_eqn++) {
		const CVEqTrm &vTerm = (*i_eqn)->getTerm();
		CVEqTrmConstIter j_term = vTerm.begin();
		(*pP)(j_term->meshID(), j_term->nodeID(), j_term->dof()) = 0.0;
		double slave_val = (*pV)(j_term->meshID(), j_term->nodeID(), j_term->dof());
		j_term++;
		for (; j_term != vTerm.end(); j_term++) {
			(*pP)(j_term->meshID(), j_term->nodeID(), j_term->dof()) +=
				j_term->coef() * slave_val;
		}
	}
}
}
