/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeFunctionCatalog.cpp
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
#include "ShapeFunctionCatalog.h"
using namespace pmw;
CShapeFunctionCatalog::CShapeFunctionCatalog()
{
    uiint itype;
    mvShapeNum.resize(ShapeType::Limit);
    for(itype=0; itype< ShapeType::Limit; itype++) {
        switch(itype) {
        case(ShapeType::Hexa81):
        case(ShapeType::Hexa82):
            mvShapeNum[itype]= 8;
            break;
        case(ShapeType::Hexa201):
        case(ShapeType::Hexa202):
        case(ShapeType::Hexa203):
            mvShapeNum[itype]= 20;
            break;
        case(ShapeType::HexaNic111):
        case(ShapeType::HexaNic118):
        case(ShapeType::HexaNic1127):
            mvShapeNum[itype]= 11;
            break;
        case(ShapeType::Tetra41):
            mvShapeNum[itype]= 4;
            break;
        case(ShapeType::Tetra101):
        case(ShapeType::Tetra104):
        case(ShapeType::Tetra1015):
            mvShapeNum[itype]= 10;
            break;
        case(ShapeType::Prism62):
            mvShapeNum[itype]= 6;
            break;
        case(ShapeType::Prism156):
        case(ShapeType::Prism159):
        case(ShapeType::Prism1518):
            mvShapeNum[itype]= 15;
            break;
        case(ShapeType::Quad41):
            mvShapeNum[itype]= 4;
            break;
        case(ShapeType::Quad84):
        case(ShapeType::Quad89):
            mvShapeNum[itype]= 8;
            break;
        case(ShapeType::Triangle31):
            mvShapeNum[itype]= 3;
            break;
        case(ShapeType::Triangle63):
            mvShapeNum[itype]= 6;
            break;
        case(ShapeType::Line21):
            mvShapeNum[itype]= 2;
            break;
        case(ShapeType::Line32):
            mvShapeNum[itype]= 3;
            break;
        }
    };
    mvIntegPointNum.resize(ShapeType::Limit);
    for(itype=0; itype< ShapeType::Limit; itype++) {
        switch(itype) {
        case(ShapeType::Hexa81):
        case(ShapeType::Hexa201):
        case(ShapeType::HexaNic111):
        case(ShapeType::Tetra41):
        case(ShapeType::Tetra101):
        case(ShapeType::Quad41):
        case(ShapeType::Triangle31):
        case(ShapeType::Line21):
            mvIntegPointNum[itype]= 1;
            break;
        case(ShapeType::Hexa82):
        case(ShapeType::Hexa202):
        case(ShapeType::HexaNic118):
            mvIntegPointNum[itype]= 8;
            break;
        case(ShapeType::Hexa203):
        case(ShapeType::HexaNic1127):
            mvIntegPointNum[itype]= 27;
            break;
        case(ShapeType::Tetra104):
        case(ShapeType::Quad84):
            mvIntegPointNum[itype]= 4;
            break;
        case(ShapeType::Tetra1015):
            mvIntegPointNum[itype]= 15;
            break;
        case(ShapeType::Prism62):
        case(ShapeType::Line32):
            mvIntegPointNum[itype]= 2;
            break;
        case(ShapeType::Prism156):
            mvIntegPointNum[itype]= 6;
            break;
        case(ShapeType::Prism159):
        case(ShapeType::Quad89):
            mvIntegPointNum[itype]= 9;
            break;
        case(ShapeType::Prism1518):
            mvIntegPointNum[itype]= 18;
            break;
        case(ShapeType::Triangle63):
            mvIntegPointNum[itype]= 3;
            break;
        default:
            break;
        }
    };
    mvIntegPair.resize(ShapeType::Limit);
    for(itype=0; itype< ShapeType::Limit; itype++) {
        mvIntegPair[itype].first= mvIntegPointNum[itype];
        mvIntegPair[itype].second= mvShapeNum[itype];
    };
    mvIntegPointNumOfElemType.resize(ElementType::Limit);
    for(itype=0; itype< ElementType::Limit; itype++) {
        switch(itype) {
        case(ElementType::Hexa):
            mvIntegPointNumOfElemType[itype].resize(2);
            mvIntegPointNumOfElemType[itype][0]= 1;
            mvIntegPointNumOfElemType[itype][1]= 8;
            break;
        case(ElementType::Hexa2):
            mvIntegPointNumOfElemType[itype].resize(3);
            mvIntegPointNumOfElemType[itype][0]=  1;
            mvIntegPointNumOfElemType[itype][1]=  8;
            mvIntegPointNumOfElemType[itype][2]= 27;
            break;
        case(ElementType::Tetra):
            mvIntegPointNumOfElemType[itype].resize(1);
            mvIntegPointNumOfElemType[itype][0]= 1;
            break;
        case(ElementType::Tetra2):
            mvIntegPointNumOfElemType[itype].resize(3);
            mvIntegPointNumOfElemType[itype][0]=  1;
            mvIntegPointNumOfElemType[itype][1]=  4;
            mvIntegPointNumOfElemType[itype][2]= 15;
            break;
        case(ElementType::Prism):
            mvIntegPointNumOfElemType[itype].resize(1);
            mvIntegPointNumOfElemType[itype][0]= 2;
            break;
        case(ElementType::Prism2):
            mvIntegPointNumOfElemType[itype].resize(3);
            mvIntegPointNumOfElemType[itype][0]=  6;
            mvIntegPointNumOfElemType[itype][1]=  9;
            mvIntegPointNumOfElemType[itype][2]= 18;
            break;
        case(ElementType::Quad):
            mvIntegPointNumOfElemType[itype].resize(1);
            mvIntegPointNumOfElemType[itype][0]= 1;
            break;
        case(ElementType::Quad2):
            mvIntegPointNumOfElemType[itype].resize(2);
            mvIntegPointNumOfElemType[itype][0]=  4;
            mvIntegPointNumOfElemType[itype][1]=  9;
            break;
        case(ElementType::Triangle):
            mvIntegPointNumOfElemType[itype].resize(1);
            mvIntegPointNumOfElemType[itype][0]= 1;
            break;
        case(ElementType::Triangle2):
            mvIntegPointNumOfElemType[itype].resize(1);
            mvIntegPointNumOfElemType[itype][0]= 3;
            break;
        case(ElementType::Beam):
            mvIntegPointNumOfElemType[itype].resize(1);
            mvIntegPointNumOfElemType[itype][0]= 1;
            break;
        case(ElementType::Beam2):
            mvIntegPointNumOfElemType[itype].resize(1);
            mvIntegPointNumOfElemType[itype][0]= 2;
            break;
        }
    };
}
CShapeFunctionCatalog::~CShapeFunctionCatalog()
{
    ;
}
integPair& CShapeFunctionCatalog::getIntegPair(const uiint& shapeType)
{
    return mvIntegPair[shapeType];
}
uiint& CShapeFunctionCatalog::NumOfShape(const uiint& shapeType)
{
    return mvShapeNum[shapeType];
}
uiint& CShapeFunctionCatalog::NumOfIntegPoint(const uiint& shapeType)
{
    return mvIntegPointNum[shapeType];
}
uiint CShapeFunctionCatalog::NumOfIntegPointType(const uiint& elemType)
{
    return mvIntegPointNumOfElemType[elemType].size();
}
uiint& CShapeFunctionCatalog::NumOfIntegPoint(const uiint& elemType, const uiint& integPtIndex)
{
    return mvIntegPointNumOfElemType[elemType][integPtIndex];
}
