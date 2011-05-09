/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ElementType.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef ELEMENT_TYPE_F9D2500C_A695_42ec_AE76_94578C5155EA
#define ELEMENT_TYPE_F9D2500C_A695_42ec_AE76_94578C5155EA
namespace pmw{
struct ElementType{
    enum{
        Polygon,
        Hexa,
        Hexa2,
        HexaNic,
        Prism,
        Prism2,
        Pyramid,
        Pyramid2,
        Tetra,
        Tetra2,
        Quad,
        Quad2,
        Triangle,
        Triangle2,
        Beam,
        Beam2,
        Point,
        Limit
    };
};
struct ShapeType{
    enum{
        Hexa81,
        Hexa82,
        Hexa201,
        Hexa202,
        Hexa203,
        HexaNic111,
        HexaNic118,
        HexaNic1127,
        Tetra41,
        Tetra101,
        Tetra104,
        Tetra1015,
        Prism62,
        Prism156,
        Prism159,
        Prism1518,
        Quad41,
        Quad84,
        Quad89,
        Triangle31,
        Triangle63,
        Line21,
        Line32,
        Limit
    };
};
struct ElementPropType{
    enum{
        Face,
        Edge,
        Node,
        Limit
    };
};
struct MaterialElement{
    enum{
        HeatTransfer,
        Spring,
        Damper,
        Limit
    };
};
struct BaseElementType{
    enum{
        Solid,
        Shell,
        Beam,
        Point,
        Cell,
        MaterialElement,
        Limit
    };
};
}
#endif 
