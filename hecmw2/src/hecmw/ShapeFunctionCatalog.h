/* 
 * File:   ShapeFunctionCatalog.h
 * Author: ktakeda
 *
 *  積分点数,節点数 の一覧
 *
 * Created on 2010/02/12, 16:26
 */
#include "TypeDef.h"
#include <utility> //pair

#include "ElementType.h"

typedef pair<uiint,uiint> integPair;

namespace pmw{
#ifndef _SHAPEFUNCTIONCATALOG_H
#define	_SHAPEFUNCTIONCATALOG_H
class CShapeFunctionCatalog{
private:
    CShapeFunctionCatalog();
public:
    static CShapeFunctionCatalog* Instance(){
        static CShapeFunctionCatalog moCatalog;
        return &moCatalog;
    }
    ~CShapeFunctionCatalog();

private:
    // {外部が必要とする情報}
    // 1.形状関数の種類数
    //  ->  struct ShapeType <= enumの終端は"Limit"
    //
    // 2.形状関数種類別-形状関数の数
    //
    // 3.形状関数種類別-積分点の数
    // 
    vuint mvShapeNum;     //形状関数種類別-形状関数の数(節点の数)
    vuint mvIntegPointNum;//形状関数種類別-積分点の数
    vector<integPair> mvIntegPair;//first=積分点数, second=節点数

    vvuint mvIntegPointNumOfElemType;//要素別-積分点種類-積分点数
public:
    // --
    // 返り値: pair.first=積分点数, pair.second=節点数
    // --
    integPair& getIntegPair(const uiint& shapeType);

    uiint& NumOfIntegPoint(const uiint& shapeType);
    uiint& NumOfShape(const uiint& shapeType);

    uiint  NumOfIntegPointType(const uiint& elemType);
    uiint& NumOfIntegPoint(const uiint& elemType, const uiint& integPtIndex);
};
#endif	/* _SHAPEFUNCTIONCATALOG_H */
}












