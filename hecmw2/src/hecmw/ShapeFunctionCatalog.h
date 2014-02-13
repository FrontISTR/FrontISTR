/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeFunctionCatalog.h
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
#include <utility>
#include "ElementType.h"
typedef pair<uiint,uiint> integPair;
namespace pmw
{
#ifndef _SHAPEFUNCTIONCATALOG_H
#define	_SHAPEFUNCTIONCATALOG_H
class CShapeFunctionCatalog
{
private:
    CShapeFunctionCatalog();
public:
    static CShapeFunctionCatalog* Instance() {
        static CShapeFunctionCatalog moCatalog;
        return &moCatalog;
    }
    ~CShapeFunctionCatalog();
private:
    vuint mvShapeNum;
    vuint mvIntegPointNum;
    vector<integPair> mvIntegPair;
    vvuint mvIntegPointNumOfElemType;
public:
    integPair& getIntegPair(const uiint& shapeType);
    uiint& NumOfIntegPoint(const uiint& shapeType);
    uiint& NumOfShape(const uiint& shapeType);
    uiint  NumOfIntegPointType(const uiint& elemType);
    uiint& NumOfIntegPoint(const uiint& elemType, const uiint& integPtIndex);
};
#endif	/* _SHAPEFUNCTIONCATALOG_H */
}
