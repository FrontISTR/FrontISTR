/* 
 * File:   CommHexa.h
 * Author: ktakeda
 *
 * Created on 2009/09/01, 15:56
 */
#include "CommElement.h"

namespace pmw{
#ifndef _COMMHEXA_H
#define	_COMMHEXA_H
class CCommHexa:public CCommElement{
public:
    CCommHexa();
    virtual ~CCommHexa();

public:
    //debug method
    virtual bool isTypeCoincidence();

    virtual uint getShapeType(){ return ElementType::Hexa;}
    virtual uint getBaseShapeType(){ return BaseElementType::Solid;}


    // prolongation
    //virtual void setupProgNodeRank(const uint& mgLevel);//ProgElemのNodeRankの決定.<= Edge,Face,Volumeのランクを決定と同義
    
};
#endif	/* _COMMHEXA_H */
}



