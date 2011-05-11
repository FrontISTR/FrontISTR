
#include "ContactNode.h"

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
    mvMasterNode.clear();
    mvSlaveNode.clear();

    mvMasterNode.resize(0);
    mvSlaveNode.resize(0);
}
COctreeKnot::~COctreeKnot()
{
    //----
    //mvChildKnotは,ContactMeshで廃棄
    //----
}



// 枝分かれする子Knotの生成
//
void COctreeKnot::createChildKnot()
{
    COctreeKnot *pKnot;
    mvChildKnot.reserve(8);
    // 枝分かれ先の生成
    uiint i;
    for(i=0; i< 8; i++){
        pKnot= new COctreeKnot;
        pKnot->setParentKnot(this);//子供に親を教えておく.
        mvChildKnot.push_back(pKnot);
    };

    double inX,inY,inZ;//中間値
    inX=(maxX + minX)*0.5;
    inY=(maxY + minY)*0.5;
    inZ=(maxZ + minZ)*0.5;
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

// 自身のItem(マスター面の頂点,スレーブ点)をChildKnotへ分配
//
void COctreeKnot::distItem()
{
    double inX,inY,inZ;//中間値
    inX=(maxX + minX)*0.5;
    inY=(maxY + minY)*0.5;
    inZ=(maxZ + minZ)*0.5;

    bool bX,bY,bZ;         //標準判定のブール
    bool beX,beY,beZ;      //等号のケースのブール
    CContactNode* pConNode;//選別するContactNode
    
    //スレーブ点の分配
    uiint numOfItem= mvSlaveNode.size();
    uiint i;
    for(i=0; i< numOfItem; i++){
        pConNode= mvSlaveNode[i];

        bX=false; bY=false; bZ=false; //標準判定
        beX=false;beY=false;beZ=false;//等号のケース

        uiint distID;//debug用途

        // X軸分割 判定
        if(pConNode->getX() < inX) bX=true;//Xは小
        // Y軸分割 判定
        if(pConNode->getY() < inY) bY=true;//Yは小
        // Z軸分割 判定
        if(pConNode->getZ() < inZ) bZ=true;//Zは小

        // ConNodeの分配(下4個)
        if( bX &&  bY &&  bZ){ mvChildKnot[0]->addSlaveNode(pConNode); distID=0;}//distIDはデバッグ用途 -> 現在,未使用
        if(!bX &&  bY &&  bZ){ mvChildKnot[1]->addSlaveNode(pConNode); distID=1;}
        if(!bX && !bY &&  bZ){ mvChildKnot[2]->addSlaveNode(pConNode); distID=2;}
        if( bX && !bY &&  bZ){ mvChildKnot[3]->addSlaveNode(pConNode); distID=3;}
        // ConNodeの分配(上4個)
        if( bX &&  bY && !bZ){ mvChildKnot[4]->addSlaveNode(pConNode); distID=4;}
        if(!bX &&  bY && !bZ){ mvChildKnot[5]->addSlaveNode(pConNode); distID=5;}
        if(!bX && !bY && !bZ){ mvChildKnot[6]->addSlaveNode(pConNode); distID=6;}
        if( bX && !bY && !bZ){ mvChildKnot[7]->addSlaveNode(pConNode); distID=7;}


        //    // 分割線上に載ってしまったConNodeを重複して配分(スレーブ点だけ)
        //    // ----
        //    if(pConNode->getX() == inX) beX=true;
        //    if(pConNode->getY() == inY) beY=true;
        //    if(pConNode->getZ() == inZ) beZ=true;
        //    // X軸分割線上
        //    if( beX ){
        //        mvChildKnot[0]->addSlaveNode(pConNode);
        //        mvChildKnot[3]->addSlaveNode(pConNode);
        //        mvChildKnot[4]->addSlaveNode(pConNode);
        //        mvChildKnot[7]->addSlaveNode(pConNode);
        //    }
        //    // Y軸分割線上
        //    if( beY ){
        //        mvChildKnot[0]->addSlaveNode(pConNode);
        //        mvChildKnot[1]->addSlaveNode(pConNode);
        //        mvChildKnot[4]->addSlaveNode(pConNode);
        //        mvChildKnot[5]->addSlaveNode(pConNode);
        //    }
        //    // Z軸分割線上
        //    if( beZ ){
        //        mvChildKnot[0]->addSlaveNode(pConNode);
        //        mvChildKnot[1]->addSlaveNode(pConNode);
        //        mvChildKnot[2]->addSlaveNode(pConNode);
        //        mvChildKnot[3]->addSlaveNode(pConNode);
        //    }
        //    // 中心点
        //    if(beX && beY && beZ){
        //        uint ii;
        //        for(ii=0; ii< 8; ii++) mvChildKnot[ii]->addSlaveNode(pConNode);
        //    }
        
    };//SlaveNodeループ


    // マスター面の頂点の分配
    numOfItem= mvMasterNode.size();
    for(i=0; i< numOfItem; i++){
        pConNode= mvMasterNode[i];

        bX=false; bY=false; bZ=false; //標準判定

        uiint distID;//debug用途

        // X軸分割 判定
        if(pConNode->getX() < inX) bX=true;//Xは小
        // Y軸分割 判定
        if(pConNode->getY() < inY) bY=true;//Yは小
        // Z軸分割 判定
        if(pConNode->getZ() < inZ) bZ=true;//Zは小

        // ConNodeの分配(下4個)
        if( bX &&  bY &&  bZ){ mvChildKnot[0]->addMasterNode(pConNode); distID=0;}//distIDはデバッグ用途 -> 現在,未使用
        if(!bX &&  bY &&  bZ){ mvChildKnot[1]->addMasterNode(pConNode); distID=1;}
        if(!bX && !bY &&  bZ){ mvChildKnot[2]->addMasterNode(pConNode); distID=2;}
        if( bX && !bY &&  bZ){ mvChildKnot[3]->addMasterNode(pConNode); distID=3;}
        // ConNodeの分配(上4個)
        if( bX &&  bY && !bZ){ mvChildKnot[4]->addMasterNode(pConNode); distID=4;}
        if(!bX &&  bY && !bZ){ mvChildKnot[5]->addMasterNode(pConNode); distID=5;}
        if(!bX && !bY && !bZ){ mvChildKnot[6]->addMasterNode(pConNode); distID=6;}
        if( bX && !bY && !bZ){ mvChildKnot[7]->addMasterNode(pConNode); distID=7;}

    };//MasterNodeループ

    //    //debug
    //    for(i=0; i < 8; i++){
    //        cout << "octreeknot distItem:: ChildKnot[" << i << "] numOfMasterNode= " << mvChildKnot[i]->getNumOfMasterNode() << endl;
    //    };
}


// ItemであるMasterNodeとSlaveNodeに
//          OctreeKnotのレイヤーとIDをセット
//
void COctreeKnot::setItemProp()
{
    CContactNode *pConNode;
    uiint numOfItem;
    uiint inode;
    
    numOfItem= mvMasterNode.size();
    for(inode=0; inode < numOfItem; inode++){
        pConNode= mvMasterNode[inode];
        pConNode->setOctreeID(mLayer, mID);
    };
    
    numOfItem= mvSlaveNode.size();
    for(inode=0; inode < numOfItem; inode++){
        pConNode= mvSlaveNode[inode];
        pConNode->setOctreeID(mLayer, mID);
    };
}






