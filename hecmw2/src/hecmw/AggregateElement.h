/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/AggregateElement.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
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
namespace pmw
{
class CAggregateElement
{
public:
    CAggregateElement();
    virtual ~CAggregateElement();
protected:
    uiint mnCoreNodeIndex;
    vector<CElement*> mvElement;
public:
    void setCoreNodeIndex(const uiint& index) {
        mnCoreNodeIndex = index;
    }
    void reserve(const uiint& res_size) {
        mvElement.reserve(res_size);
    }
    void push(CElement* pElement) {
        mvElement.push_back(pElement);
    }
    CElement* get(const uiint& index) {
        return mvElement[index];
    }
    uiint getNumOfElement() {
        return mvElement.size();
    }
};
}
#endif	/* _AGGREGATEELEMENT_H */
