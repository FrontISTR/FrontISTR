/* 
 * File:   BoundaryGroup.h
 * Author: ktakeda
 *
 * Boundary Template
 *
 *
 * Created on 2009/05/21, 14:17
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
    map<unsigned int, unsigned int, less<unsigned int> > mIndexMap;// Meshの要素IndexからmvBoundaryのIndex番号を取得

public:
    // BoundaryMeshオブジェクトの確保
    // --
    void reserve(const unsigned int& res_size){ mvBoundary.reserve(res_size); }
    void push(T pBoundary){ 
        
        mvBoundary.push_back(pBoundary);
        
        unsigned int index = mvBoundary.size() - 1;
        unsigned int id = pBoundary->getID();
        
        mIndexMap[id] = index;
    }


    // 各Boundaryオブジェクトを生成した際にIDはセットするのでコメントアウト.
    //
    //void setID(const uint& index, const uint& id){
    //    mvBoundary[index]->setID(id);
    //}

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
    // object 数
    unsigned int NumOfBoundary(){
        return mvBoundary.size();
    }
};
}
#endif	/* _BOUNDARYGROUP_H */

