/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ElementType.h
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
#ifndef ELEMENT_TYPE_F9D2500C_A695_42ec_AE76_94578C5155EA
#define ELEMENT_TYPE_F9D2500C_A695_42ec_AE76_94578C5155EA
namespace pmw
{
struct ElementType {
    enum {
        Polygon=0,
        Hexa=1,
        Hexa2=2,
        HexaNic=3,
        Prism=4,
        Prism2=5,
        Tetra=6,
        Tetra2=7,
        Quad=8,
        Quad2=9,
        Triangle=10,
        Triangle2=11,
        Beam=12,
        Beam2=13,
        Line=14,
        Line2=15,
        Point=16,
        Limit=17
    };
};
struct FistrElementType {
    enum {
        Polygon=1,
        Point,
        Beam=111,       //ロッド、リンク、トラス(1次)
        Beam2=112,      //ロッド、リンク、トラス(2次)
        Triangle=231,   //三角形要素(1次)
        Triangle2=232,  //三角形要素(2次)
        Quad=241,       //四辺形要素(1次)
        Quad2=242,      //四辺形要素(2次)
        Tetra=341,      //四面体要素(1次)
        Tetra2=342,     //四面体要素(2次)
        Prism=351,      //三角柱要素(1次)
        Prism2=352,     //三角柱要素(2次)
        Hexa=361,       //六面体要素(1次)
        Hexa2=362,      //六面体要素(2次)
        IFaceQuad=541,  //インターフェース要素(四辺形断面、1次)
        IFaceQuad2=542, //インターフェース要素(四辺形断面、2次)
        TriShell=731,   //三角形シェル要素(1次)
        TriShell2=732,  //三角形シェル要素(2次)
        QuadShell=741,  //四辺形シェル要素(1次)
        QuadShell2=742, //四辺形シェル要素(2次)
        Limit
    };
};
struct ElementOrder {
    enum {
        Zero=0,
        First=1,
        Second=2,
        Limit=3
    };
};
struct ShapeType {
    enum {
        Hexa81=0,
        Hexa82=1,
        Hexa201=2,
        Hexa202=3,
        Hexa203=4,
        HexaNic111=5,
        HexaNic118=6,
        HexaNic1127=7,
        Tetra41=8,
        Tetra101=9,
        Tetra104=10,
        Tetra1015=11,
        Prism62=12,
        Prism156=13,
        Prism159=14,
        Prism1518=15,
        Quad41=16,
        Quad84=17,
        Quad89=18,
        Triangle31=19,
        Triangle63=20,
        Line21=21,
        Line32=22,
        Limit=23
    };
};
struct ElementPropType {
    enum {
        Face=0,
        Edge=1,
        Node=2,
        Limit=3
    };
};
struct MaterialElement {
    enum {
        HeatTransfer=0,
        Spring=1,
        Damper=2,
        Limit=3
    };
};
struct BaseElementType {
    enum {
        Solid=0,
        Shell=1,
        Beam=2,
        Point=3,
        Cell=4,
        MaterialElement=5,
        Limit=6
    };
};
}
#endif
