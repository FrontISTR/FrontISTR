/* 
 * File:   OctreeKnot.h
 * Author: ktakeda
 *
 * Octreeの節:List構造
 *   -> 接合面メッシュ(MPCメッシュ)の検索絞り込み用途
 *
 *
 *  [Knot-Position]
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
 *
 * Created on 2009/11/26, 19:37
 */
#include "TypeDef.h"
#include "ContactNode.h"//アイテム

#include <iostream>

namespace pmw{
#ifndef _OCTREEKNOT_H
#define	_OCTREEKNOT_H
class COctreeKnot{
public:
    COctreeKnot();
    virtual ~COctreeKnot();

protected:
    uint mLayer;//レイヤー番号
    uint mID;   //箱番号(OctreeBox番号)
    
    
    COctreeKnot* mpParentKnot;       //親のKnot
    vector<COctreeKnot*> mvChildKnot;//子供のKnot:配列Index => Position

    vector<CContactNode*> mvMasterNode;//マスター面の頂点のコンテナ
    vector<CContactNode*> mvSlaveNode; //スレーブ点のコンテナ

    // Itemの範囲(BoundingBox):3軸の範囲
    double minX,maxX;
    double minY,maxY;
    double minZ,maxZ;

public:
    //レイヤー,箱番号(OctreeBox番号)
    void setLayerID(const uint& layer){ mLayer= layer;}
    void setID(const uint& id){ mID= id;}
    uint& getLayerID(){ return mLayer;}
    uint& getID(){ return mID;}

    //アイテムにOctreeKnotのレイヤーとIDをセット
    void setItemProp();

    //// アイテム(ContactNode*) 入力アクセッサー
    //void reserveItem(const uint& res_size){ mvItem.reserve(res_size);}
    //void resizeItem(const uint& res_size){ mvItem.resize(res_size);}
    //void setItem(CContactNode* pConNode, const uint& index){ mvItem[index]= pConNode;}
    //void addItem(CContactNode* pConNode){ mvItem.push_back(pConNode);}
    //// アイテム 出力アクセッサー
    //uint getNumOfItem(){ return mvItem.size();}
    //vector<CContactNode*>& getItem(){ return mvItem;}
    //CContactNode*  getItem(const uint& index){ return mvItem[index];}

    // マスター面の頂点(ContactNode*)
    void reserveMasterNode(const uint& res_size){ mvMasterNode.reserve(res_size);}
    void addMasterNode(CContactNode* pConNode){ mvMasterNode.push_back(pConNode);}
    uint getNumOfMasterNode(){ return mvMasterNode.size();}
    vector<CContactNode*>& getMasterNode(){ return mvMasterNode;}
    CContactNode*  getMasterNode(const uint& index){ return mvMasterNode[index];}

    // スレーブ点
    void reserveSlaveNode(const uint& res_size){ mvSlaveNode.reserve(res_size);}
    void addSlaveNode(CContactNode* pConNode){ mvSlaveNode.push_back(pConNode);}
    uint getNumOfSlaveNode(){ return mvSlaveNode.size();}
    vector<CContactNode*>& getSlaveNode(){ return mvSlaveNode;}
    CContactNode*  getSlaveNode(const uint& index){ return mvSlaveNode[index];}


    // Knotの場合の範囲(Box)指定
    void setX(const double& min, const double& max){ minX=min; maxX=max;}
    void setY(const double& min, const double& max){ minY=min; maxY=max;}
    void setZ(const double& min, const double& max){ minZ=min; maxZ=max;}
    // Knot範囲の提供
    double& getMinX(){return minX;} double& getMaxX(){return maxX;}
    double& getMinY(){return minY;} double& getMaxY(){return maxY;}
    double& getMinZ(){return minZ;} double& getMaxZ(){return maxZ;}


    // 親Knot：自身(Knot)の上流
    void setParentKnot(COctreeKnot* pKnot){ mpParentKnot= pKnot;}
    COctreeKnot* getParentKnot(){ return mpParentKnot;}

    // 子Knot：枝分かれ先のKnot
    // --
    void createChildKnot();
    COctreeKnot* getChildKnot(const uint& pos){ return mvChildKnot[pos];}
    void distItem();//自身のItem(マスター面の頂点,スレーブ点)をChildKnotへ分配
};
#endif	/* _OCTREEKNOT_H */
}















