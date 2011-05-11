/* 
 * File:   AggregateElement.h
 * Author: ktakeda
 *
 * Created on 2009/05/29, 15:50
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
    uiint mnCoreNodeIndex;//AggElementの中心にいるNodeのIndex
    vector<CElement*> mvElement;

public:
    void setCoreNodeIndex(const uiint& index){ mnCoreNodeIndex = index;}//中心に存在する節点Index

    void reserve(const uiint& res_size){ mvElement.reserve(res_size);}
    void push(CElement* pElement){ mvElement.push_back(pElement);}
    CElement* get(const uiint& index){ return mvElement[index];}

    uiint getNumOfElement(){ return mvElement.size();}
};
}

#endif	/* _AGGREGATEELEMENT_H */

