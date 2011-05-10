/*
 * ContactMatrix.h
 *
 *  Created on: Oct 16, 2009
 *      Author: goto
 */

#ifndef CONTACTMATRIX_H_
#define CONTACTMATRIX_H_

namespace pmw
{

class CAssyVector;

class CContactMatrix
{
public:
	CContactMatrix();
	virtual ~CContactMatrix();
	void multVectorAdd(CAssyVector *pP, CAssyVector *pX);
};

}

#endif /* CONTACTMATRIX_H_ */
