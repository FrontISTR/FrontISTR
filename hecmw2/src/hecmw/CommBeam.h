/* 
 * File:   CommBeam.h
 * Author: ktakeda
 *
 * Created on 2009/09/10, 16:35
 */
#include "CommElement.h"

namespace pmw{
#ifndef _COMMBEAM_H
#define	_COMMBEAM_H
class CCommBeam:public CCommElement{
public:
    CCommBeam();
    virtual ~CCommBeam();
    
public:
    //debug method
    virtual bool isTypeCoincidence();

    virtual uiint getShapeType(){ return ElementType::Beam;}
    virtual uiint getBaseShapeType(){ return BaseElementType::Beam;}

    //virtual void setupProgNodeRank(const uint& mgLevel);//ProgElemのNodeRankの決定.<= Edge,Face,Volumeのランクを決定と同義
    
};
#endif	/* _COMMBEAM_H */
}

