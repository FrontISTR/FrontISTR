/*
 * MPCMatrix.h
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
 */

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
		int id = mvEquation.size();
		mnEquation = id + 1;
		mvEquation.push_back(equation);
	};

private:
	size_t mnEquation;
	std::vector<CEquation*> mvEquation;
};
#endif /* MPCMATRIX_H_ */
}


