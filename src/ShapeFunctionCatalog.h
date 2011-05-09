/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ShapeFunctionCatalog.h
|
|                     Written by T.Takeda,    2010/06/01
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
typedef pair<uint,uint> integPair;
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
    vuint mvShapeNum;     
    vuint mvIntegPointNum;
    vector<integPair> mvIntegPair;
    vvuint mvIntegPointNumOfElemType;
public:
    integPair& getIntegPair(const uint& shapeType);
    uint& NumOfIntegPoint(const uint& shapeType);
    uint& NumOfShape(const uint& shapeType);
    uint  NumOfIntegPointType(const uint& elemType);
    uint& NumOfIntegPoint(const uint& elemType, const uint& integPtIndex);
};
#endif	/* _SHAPEFUNCTIONCATALOG_H */
}
