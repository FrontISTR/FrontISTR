/* 
 * File:   BoundingBox.h
 * Author: ktakeda
 *
 * MPCマスター面周囲のスレーブ点を絞り込み検索
 *
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
 *  Z             ABB
 *  |   .Y
 *  | .
 *  o ----- X
 *
 *
 *
 *           Z
 *   ________|_________
 *  |        |         |
 *  |         ---->X   |
 *   ------------------
 *          OBB 
 *
 *
 *            2009.12.15
 * Created on 2009/10/14, 13:32
 */
#include "TypeDef.h"

#include "SkinFace.h"//接合面(スキン):マスター,スレーブ面

namespace pmw{
#ifndef _BOUNDINGBOX_H
#define	_BOUNDINGBOX_H
class CBoundingBox{
public:
    CBoundingBox();
    virtual ~CBoundingBox();

protected:
    //ABB
    vdouble mMinCoord,mMaxCoord;//ABB範囲

    //OBB
    vdouble mvOBBCenter;//OBB座標の中心点(x,y,z)
    vdouble mvOBBX;//OBB座標のX方向ベクトル
    vdouble mvOBBY;//OBB座標のY方向ベクトル
    vdouble mvOBBZ;//OBB座標のZ方向ベクトル
    vdouble mvE;//OBB座標のX,Y,Z方向の範囲
public:
    //ABB範囲
    void sizingABB(CSkinFace* pFace);//ABBの最大,最小座標：マスター面のサイズを取得して決定
    
    //ABB判定
    bool judgeABB(CSkinFace* pFace);//SkinFace(スレーブ面)の判定
    bool judgeABB(CContactNode* pConNode);//ContactNode(スレーブ面ConNode)の判定
    
    
    
    //OBB座標軸,範囲
    void sizingOBB(CSkinFace* pFace);//OBBの座標軸方向と範囲：マスター面の向きとサイズから決定

    //OBB判定
    bool judgeOBB(CContactNode* pConNode);//ContactNode(スレーブ面ConNode)の判定
};
#endif	/* _BOUNDINGBOX_H */
}



