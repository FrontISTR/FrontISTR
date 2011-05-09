/* 
 * File:   CommQuad.h
 * Author: ktakeda
 *
 * Created on 2009/09/09, 18:31
 */
#include "CommElement.h"

namespace pmw{
#ifndef _COMMQUAD_H
#define	_COMMQUAD_H
class CCommQuad:public CCommElement{
public:
    CCommQuad();
    virtual ~CCommQuad();

public:
    //debug method
    virtual bool isTypeCoincidence();

    virtual uint getShapeType(){ return ElementType::Quad;}
    virtual uint getBaseShapeType(){ return BaseElementType::Shell;}

    //virtual void setupProgNodeRank(const uint& mgLevel);//ProgElemのNodeRankの決定.<= Edge,Face,Volumeのランクを決定と同義
    
};
#endif	/* _COMMQUAD_H */
}

