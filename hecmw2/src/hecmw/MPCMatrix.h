/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/MPCMatrix.h
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
#include "TypeDef.h"
#include <vector>
namespace pmw
{
class CAssyVector;
class CEquation;
#ifndef MPCMATRIX_H_
#define MPCMATRIX_H_
class CMPCMatrix
{
public:
	CMPCMatrix();
	virtual ~CMPCMatrix();
	void multVector(CAssyVector *pV, CAssyVector *pP) const;
	void transMultVector(CAssyVector *pV, CAssyVector *pP) const;
	void addEquation(CEquation* equation){
		uiint id = mvEquation.size();
		mnEquation = id + 1;
		mvEquation.push_back(equation);
	};
        void dump();
private:
	uiint mnEquation;
	std::vector<CEquation*> mvEquation;
};
#endif /* MPCMATRIX_H_ */
}
