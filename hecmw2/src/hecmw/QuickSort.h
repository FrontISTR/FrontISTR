/* 
 * File:   QuickSort.h
 * Author: ktakeda
 *
 * Created on 2009/09/15, 19:49
 */
#include "TypeDef.h"
#include <iostream>

namespace pmw{
#ifndef _QUICKSORT_H
#define	_QUICKSORT_H
//---
// QuickSortID
//---
// --
// CNode*, CElement* を T としてIDを比較してソート
// --
template<class T>
void QuicksortID(vector<T>& val, const uiint& start, const uiint& end){
    uiint i= start;//開始インデックス
    uiint j= end;  //終了インデックス
    uiint base= (val[start]->getID() + val[end]->getID())/2; //基準値(平均値)
    //int base= val[start]->getID();//基準値

    ////debug
    //std::cout << "QuicksortID  base => " << base << std::endl;

    while(true){
        while(val[i]->getID() < base) i++; //基準値より大きい値を前から探す
        while(val[j]->getID() > base) j--; //基準値より小さい値を後ろから探す

        if(i>=j) break;//交差したらブレーク(交換すべき値がない)

        //値の交換
        T temp = val[i];
        val[i] = val[j];
        val[j] = temp;
        
        i++;//交換した一つ後ろにインデックスを進める
        j--;//交換した一つ前にインデックスを進める
    };

    //再帰処理::配列を二分して,それぞれに再帰処理
    if(start < i-1) QuicksortID(val, start, i-1);//基準値より小さいグループ
    if(end   > j+1) QuicksortID(val, j+1,   end);//基準値より大きいグループ
};
//--
// ソート ID
//--
template<class T>
void sortID(vector<T>& val, const uiint& length)
{
    uiint i;
    uiint j;

    for(i=0; i< length-1; i++){
        for(j=i+1; j< length; j++){
            //値の交換
            if(val[i]->getID() > val[j]->getID()){
                T temp = val[i];
                val[i] = val[j];
                val[j] = temp;
            }
        };
    };
    ////debug
    //for(i=0; i< length; i++)
    //    std::cout << "val[" << i << "]=> " << val[i]->getID() << std::endl;
};


// --
// ↓ 旧
// --
//template<class T>
//void QuicksortID(vector<T>& val, const uint& left, const uint& right){
//    if(left==right) return;
//
//    int  pivpos;
//    uint parpos;
//    // pivot計算 {配列の軸位置}
//    pivpos = pivotIndex(val,left,right);
//
//    //debug
//    std::cout << "QuicksortID::pivot pos => " << pivpos << std::endl;
//
//    if(pivpos != -1){
//        parpos= partitionIndex(val, left, right, val[pivpos]);//partition分割位置
//
//        //debug
//        std::cout << "QuicksortID::partition pos => " << parpos << std::endl;
//
//        QuicksortID(val, left, parpos-1);
//        QuicksortID(val, parpos, right);
//    }
//
//}
//// --
//// pivot位置の選択
//// --
//// 最初にみつかった２つの異なる要素のうち,大きいほうのIndex番号を返す.
//template<class T>
//int pivotIndex(vector<T>& val, const uint& i, const uint& j){
//
//    uint k;
//    k = i + 1;
//
//    while(k <= j && val[i]->getID()==val[k]->getID()) k++;
//    if(k > j) return -1;
//    if(val[i]->getID() >= val[k]->getID()) return i;
//
//    return k;
//};
//// --
//// partition分割位置
//// --
//template<class T>
//uint partitionIndex(vector<T>& val, const uint& i, const uint& j, T& pivotVal){
//    uint left,right;
//    left= i; right= j;
//
//    //debug
//    std::cout << "partitionIndex::pivotVal->getID()=>" << pivotVal->getID() << std::endl;
//    std::cout << "partitionIndex::" << "left=> " << left << ", right=> " << right << std::endl;
//
//    // --
//    // 検索が交差するまで繰り返し
//    while(left <= right){
//        // 軸要素以上のデータを検索
//        while(left  <= j && val[left]->getID()  <  pivotVal->getID()) left++;
//        //debug
//        std::cout << "partitionIndex::while left++ の結果 => " << left << std::endl;
//
//        // 軸要素未満のデータを検索
//        while(right >= i && val[right]->getID() >= pivotVal->getID()) right--;
//        //debug
//        std::cout << "partitionIndex::while right-- の結果 => " << right << std::endl;
//
//        if(left > right) break;
//        //入れ替え(swap)
//        T temp    = val[left];
//        val[left] = val[right];
//        val[right]= temp;
//        left++; right--;
//    };
//    return left;
//}

#endif	/* _QUICKSORT_H */
}






















