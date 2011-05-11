/* 
 * File:   CommPrism.h
 * Author: ktakeda
 *
 * Created on 2009/09/08, 14:43
 */
#include "CommElement.h"
#include "CommHexa.h"

namespace pmw{
#ifndef _COMMPRISM_H
#define	_COMMPRISM_H
class CCommPrism:public CCommElement{
public:
    CCommPrism();
    virtual ~CCommPrism();

public:
    //debug method
    virtual bool isTypeCoincidence();

    virtual uiint getShapeType(){ return ElementType::Prism;}
    virtual uiint getBaseShapeType(){ return BaseElementType::Solid;}

    //virtual void setupProgNodeRank(const uint& mgLevel);//ProgElemのNodeRankの決定.<= Edge,Face,Volumeのランクを決定と同義
    
};
#endif	/* _COMMPRISM_H */
}




