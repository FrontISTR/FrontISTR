/* 
 * File:   CommPyramid.h
 * Author: ktakeda
 *
 * Created on 2009/09/09, 18:11
 */
#include "CommElement.h"
#include "CommHexa.h"

namespace pmw{
#ifndef _COMMPYRAMID_H
#define	_COMMPYRAMID_H
class CCommPyramid:public CCommElement{
public:
    CCommPyramid();
    virtual ~CCommPyramid();

public:
    //debug method
    virtual bool isTypeCoincidence();

    virtual uint getShapeType(){ return ElementType::Pyramid;}
    virtual uint getBaseShapeType(){ return BaseElementType::Solid;}

    //virtual void setupProgNodeRank(const uint& mgLevel);//ProgElemのNodeRankの決定.<= Edge,Face,Volumeのランクを決定と同義
    
};
#endif	/* _COMMPYRAMID_H */
}

