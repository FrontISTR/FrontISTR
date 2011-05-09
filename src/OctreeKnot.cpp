//
//  OctreeKnot.cpp
//
//
//
//                  2009.11.26
//                  2009.11.26
//                  k.Takeda
#include "OctreeKnot.h"
using namespace pmw;

// construct & destruct
//
COctreeKnot::COctreeKnot()
{
    ;
}
COctreeKnot::~COctreeKnot()
{
    ;
}


// 枝分かれする子Knotの生成
//
void COctreeKnot::createChildKnot()
{
    mvChildKnot.reserve(8);
    // 枝分かれ先の生成
    uint i;
    for(i=0; i< 8; i++){
        mvChildKnot.push_back(new COctreeKnot);
    };

    double inX,inY,inZ;//中間値
    inX=(minX + maxX)*0.5; inY=(minY + maxY)*0.5; inZ=(minZ + maxZ)*0.5;
    /*
     *   枝分かれ先の範囲
     *
     *       7.........6
     *     .         . |
     *   4 ---------5  |
     *   |          |  |
     *   |    3     |  | 2
     *   |          | .
     *   |          |
     *   0 --------- 1
     *
     *  Z
     *  |   .Y
     *  | .
     *  o ----- X
     */
    // position 0
    mvChildKnot[0]->setX(minX,inX); mvChildKnot[0]->setY(minY,inY); mvChildKnot[0]->setZ(minZ,inZ);
    // position 1
    mvChildKnot[1]->setX(inX,maxX); mvChildKnot[1]->setY(minY,inY); mvChildKnot[1]->setZ(minZ,inZ);
    // position 2
    mvChildKnot[2]->setX(inX,maxX); mvChildKnot[2]->setY(inY,maxY); mvChildKnot[2]->setZ(minZ,inZ);
    // position 3
    mvChildKnot[3]->setX(minX,inX); mvChildKnot[3]->setY(inY,maxY); mvChildKnot[3]->setZ(minZ,inZ);

    // position 4
    mvChildKnot[4]->setX(minX,inX); mvChildKnot[4]->setY(minY,inY); mvChildKnot[4]->setZ(inZ,maxZ);
    // position 5
    mvChildKnot[5]->setX(inX,maxX); mvChildKnot[5]->setY(minY,inY); mvChildKnot[5]->setZ(inZ,maxZ);
    // position 6
    mvChildKnot[6]->setX(inX,maxX); mvChildKnot[6]->setY(inY,maxY); mvChildKnot[6]->setZ(inZ,maxZ);
    // position 7
    mvChildKnot[7]->setX(minX,inX); mvChildKnot[7]->setY(inY,maxY); mvChildKnot[7]->setZ(inZ,maxZ);
}

// 自身のItemをChildKnotへ分配
//
void COctreeKnot::distItem()
{
    double inX,inY,inZ;//中間値
    inX=(minX + maxX)*0.5; inY=(minY + maxY)*0.5; inZ=(minZ + maxZ)*0.5;

    bool bX,bY,bZ;
    CContactNode* pConNode;//選別するContactNode
    uint numOfItem= mvItem.size();
    uint i;
    for(i=0; i< numOfItem; i++){
        pConNode= mvItem[i];

        bX=false; bY=false; bZ=false;

        // X軸分割 判定
        if(pConNode->getX() < inX) bX=true;//Xは小
        // Y軸分割 判定
        if(pConNode->getY() < inY) bY=true;//Yは小
        // Z軸分割 判定
        if(pConNode->getZ() < inZ) bZ=true;//Zは小

        // ConNodeの分配(下4個)
        if(bX &&  bY && bZ) mvChildKnot[0]->addItem(pConNode);
        if(!bX &&  bY && bZ) mvChildKnot[1]->addItem(pConNode);
        if(!bX && !bY && bZ) mvChildKnot[2]->addItem(pConNode);
        if(bX && !bY && bZ) mvChildKnot[3]->addItem(pConNode);
        // ConNodeの分配(上4個)
        if(bX &&  bY && !bZ) mvChildKnot[4]->addItem(pConNode);
        if(!bX &&  bY && !bZ) mvChildKnot[5]->addItem(pConNode);
        if(!bX && !bY && !bZ) mvChildKnot[6]->addItem(pConNode);
        if(bX && !bY && !bZ) mvChildKnot[7]->addItem(pConNode);
        
    };//Itemループ
}









