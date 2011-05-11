/* 
 * File:   CommTetra.h
 * Author: ktakeda
 *
 * Created on 2009/09/08, 14:13
 */
#include "CommElement.h"
#include "CommHexa.h"

namespace pmw{
#ifndef _COMMTETRA_H
#define	_COMMTETRA_H
class CCommTetra:public CCommElement{
public:
    CCommTetra();
    virtual ~CCommTetra();

public:
    //debug method
    virtual bool isTypeCoincidence();

    virtual uiint getShapeType(){ return ElementType::Tetra;}
    virtual uiint getBaseShapeType(){ return BaseElementType::Solid;}

    //virtual void setupProgNodeRank(const uint& mgLevel);//ProgElemのNodeRankの決定.<= Edge,Face,Volumeのランクを決定と同義
    
};
#endif	/* _COMMTETRA_H */
}

