/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryGroup.h
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
#include "TypeDef.h"
#include <vector>
#include <map>
using namespace std;
namespace pmw
{
#ifndef _BOUNDARYGROUP_H_5a9ad540_ba85_4b69_96c1_24b7f9f6ad67
#define	_BOUNDARYGROUP_H_5a9ad540_ba85_4b69_96c1_24b7f9f6ad67
template<class T>
class BoundaryGroup
{
protected:
    vector<T> mvBoundary;
    map<uiint, uiint, less<uiint> > mIndexMap;
public:
    void reserve(const uiint& res_size) {
        mvBoundary.reserve(res_size);
    }
    void push(T pBoundary) {
        mvBoundary.push_back(pBoundary);
        uiint index = mvBoundary.size() - 1;
        uiint id = pBoundary->getID();
        mIndexMap[id] = index;
    }
    T get_withID(const uiint& id) {
        uiint index;
        index = mIndexMap[id];
        return mvBoundary[index];
    }
    T get_withIndex(const uiint& index) {
        return mvBoundary[index];
    }
    void setupIndexMap() {
        uiint i, id;
        for(i=0; i< mvBoundary.size(); i++) {
            id = mvBoundary[i]->getID();
            mIndexMap[id] = i;
        };
    }
    uiint NumOfBoundary() {
        return mvBoundary.size();
    }
    void clear() {
        uiint nNumOfBoundary= mvBoundary.size();
        for(uiint i=0; i < nNumOfBoundary; i++) {
            delete mvBoundary[i];
        }
        mvBoundary.clear();
        mIndexMap.clear();
    }
};
#endif	/* _BOUNDARYGROUP_H */
}

