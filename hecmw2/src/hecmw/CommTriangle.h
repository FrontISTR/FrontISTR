/* 
 * File:   CommTriangle.h
 * Author: ktakeda
 *
 * Created on 2009/09/10, 16:28
 */
#include "CommElement.h"

namespace pmw{
#ifndef _COMMTRIANGLE_H
#define	_COMMTRIANGLE_H
class CCommTriangle:public CCommElement{
public:
    CCommTriangle();
    virtual ~CCommTriangle();

public:
    //debug method
    virtual bool isTypeCoincidence();

    virtual uint getShapeType(){ return ElementType::Triangle;}
    virtual uint getBaseShapeType(){ return BaseElementType::Shell;}

    //virtual void setupProgNodeRank(const uint& mgLevel);//ProgElemのNodeRankの決定.<= Edge,Face,Volumeのランクを決定と同義
    
};
#endif	/* _COMMTRIANGLE_H */
}

