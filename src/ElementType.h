//
//  Element Type
//
//  終端として,Limitを追加 2010.02.12
//
//				2010.02.12
//				2008.12.04
//				k.Takeda
#ifndef ELEMENT_TYPE_F9D2500C_A695_42ec_AE76_94578C5155EA
#define ELEMENT_TYPE_F9D2500C_A695_42ec_AE76_94578C5155EA


//#define HEXA_TYPE     0
//#define PRISM_TYPE    1
//#define TETRA_TYPE    2
//#define QUAD_TYPE     3
//#define TRIANGLE_TYPE 4
//#define BEAM_TYPE     5
//#define BAR_TYPE      6
//#define ROD_TYPE      7
//#define MASS_TYPE     8

namespace pmw{
// Element
struct ElementType{
    enum{
        Polygon,//Film
        Hexa,
        Hexa2,
        HexaNic,
        Prism,
        Prism2,
//        Pyramid,
//        Pyramid2,
        Tetra,
        Tetra2,
        Quad,
        Quad2,
        Triangle,
        Triangle2,
        Beam,
        Beam2,
        //Bar,
        //Rod,
        //Mass,
        Point,//SkinFaceの型として利用
        Limit//enumの終端
    };
};

// FrontISTR要素タイプ
//
struct FistrElementType{
    enum{
        Polygon=1,//Film
        Hexa=361,
        Hexa2=362,
        Prism=351,
        Prism2=352,
        Tetra=341,
        Tetra2=342,
        Quad=241,
        Quad2=242,
        Triangle=231,
        Triangle2=232,
        Beam,
        Beam2,
        //Bar,
        //Rod,
        //Mass,
        Point,//SkinFaceの型として利用
        Limit//enumの終端
    };
};

struct ElementOrder{
    enum{
        Zero,
        First,
        Second,
        Limit//enumの終端
    };
};

// ShapeFunctionType
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


// Property
struct ElementPropType{
    enum{
        Face,
        Edge,
        Node,
        Limit//enumの終端
    };
};

// Material
struct MaterialElement{
    enum{
        HeatTransfer,
        Spring,
        Damper,
        Limit//enumの終端
    };
};

// Base Type
struct BaseElementType{
    enum{
        Solid,
        Shell,
        Beam,
        Point,
        Cell,
        MaterialElement,
        Limit//enumの終端
    };
};
}
#endif // ELEMENT_TYPE_F9D2500C_A695_42ec_AE76_94578C5155EA


