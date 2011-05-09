/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BoundaryGroup.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _BOUNDARYGROUP_H_5a9ad540_ba85_4b69_96c1_24b7f9f6ad67
#define	_BOUNDARYGROUP_H_5a9ad540_ba85_4b69_96c1_24b7f9f6ad67
#include <vector>
#include <map>
using namespace std;
namespace pmw{
template<class T>
class BoundaryGroup{
protected:
    vector<T> mvBoundary;
    map<unsigned int, unsigned int, less<unsigned int> > mIndexMap;
public:
    void reserve(const unsigned int& res_size){ mvBoundary.reserve(res_size); }
    void push(T pBoundary){ mvBoundary.push_back(pBoundary);}
    T get_withID(const unsigned int& id){
        unsigned int index;
        index = mIndexMap[id];
        return mvBoundary[index];
    }
    T get_withIndex(const unsigned int& index){
        return mvBoundary[index];
    }
    void setupIndexMap(){
        unsigned int i, id;
        for(i=0; i< mvBoundary.size(); i++){
            id = mvBoundary[i]->getID();
            mIndexMap[id] = i;
        };
    }
    unsigned int NumOfBoundary(){
        return mvBoundary.size();
    }
};
}
#endif	/* _BOUNDARYGROUP_H */
