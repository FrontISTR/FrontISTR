/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Equation.h
|
|                     Written by T.Takeda,    2010/06/01
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
	void setSlave(int meshID, int nodeID, int dof)
	{
		clear();
		mnTerm = 1;
		mvTerm.resize(mnTerm);
		mvTerm[0].set(meshID, nodeID, dof, 1.0);
	}
	void addMaster(int meshID, int nodeID, int dof, double coef)
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
	const CEqnTerm &getTerm(size_t idx) const { return mvTerm[idx]; }
	const std::vector<CEqnTerm> &getTerm() const { return mvTerm; }
	void setConstTerm(double const_term) { mConstTerm = const_term; }
	const double &getConstTerm() const { return mConstTerm; }
	const size_t &numTerm() const { return mnTerm; }
	void clear()
	{
		mnTerm = 0;
		mvTerm.resize(0);
		mConstTerm = 0.0;
	}
private:
	size_t mnTerm;
	std::vector<CEqnTerm> mvTerm;
	double mConstTerm;
};
}
#endif /* EQUATION_H_ */
