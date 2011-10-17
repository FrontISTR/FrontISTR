/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Equation.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef EQUATION_H_
#define EQUATION_H_
#include <vector>
#include <iostream>
#include "EqnTerm.h"
namespace pmw
{
class CEquation
{
public:
	CEquation() : mnTerm(0), mConstTerm(0.0) {}
	virtual ~CEquation() {}
	void setSlave(uiint meshID, uiint nodeID, uiint dof)
	{
		clear();
		mnTerm = 1;
		mvTerm.resize(mnTerm);
		mvTerm[0].set(meshID, nodeID, dof, 1.0);
	}
	void addMaster(uiint meshID, uiint nodeID, uiint dof, double coef)
	{
		if (mnTerm < 1) {
			std::cerr << "ERROR: CEquation::addMaster: call setSlave first" << std::endl;
			exit(1);
		}
		int id = mvTerm.size();
		mnTerm = id + 1;
		mvTerm.resize(mnTerm);
		mvTerm[id].set(meshID, nodeID, dof, coef);
	}
	const CEqnTerm &getTerm(uiint idx) const { return mvTerm[idx]; }
	const std::vector<CEqnTerm> &getTerm() const { return mvTerm; }
	void setConstTerm(double const_term) { mConstTerm = const_term; }
	const double &getConstTerm() const { return mConstTerm; }
	const uiint &numTerm() const { return mnTerm; }
	void clear()
	{
		mnTerm = 0;
		mvTerm.resize(0);
		mConstTerm = 0.0;
	}
private:
	uiint mnTerm;
	std::vector<CEqnTerm> mvTerm;
	double mConstTerm;
};
}
#endif /* EQUATION_H_ */
