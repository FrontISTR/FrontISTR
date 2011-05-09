/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   AggregateElement.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _AGGREGATEELEMENT_H_7ad7e855_66ab_4e03_a3ae_d43f6d24a78d
#define	_AGGREGATEELEMENT_H_7ad7e855_66ab_4e03_a3ae_d43f6d24a78d
#include "Element.h"
namespace pmw{
class CAggregateElement{
public:
    CAggregateElement();
    virtual ~CAggregateElement();
protected:
    uint mnCoreNodeIndex;
    vector<CElement*> mvElement;
public:
    void setCoreNodeIndex(const uint& index){ mnCoreNodeIndex = index;}
    void reserve(const uint& res_size){ mvElement.reserve(res_size);}
    void push(CElement* pElement){ mvElement.push_back(pElement);}
    CElement* get(const uint& index){ return mvElement[index];}
    uint getNumOfElement(){ return mvElement.size();}
};
}
#endif	/* _AGGREGATEELEMENT_H */
