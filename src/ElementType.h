//
//  Element Type
//
//				2009.05.11
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
union ElementType{
    enum{
        Polyhedral,
        Hexa,
        Hexa2,
        Prism,
        Prism2,
        Pyramid,
        Pyramid2,
        Tetra,
        Tetra2,
        Polygon,
        Quad,
        Quad2,
        Triangle,
        Triangle2,
        Beam,
        Beam2,
        Bar,
        Rod,
        Mass
    };
};

// Property
union ElementPropType{
    enum{
        Face,
        Edge,
        Node
    };
};

// Material
union MaterialElement{
    enum{
        HeatTransfer,
        Spring,
        Damper
    };
};

// Base Type
union BaseElementType{
    enum{
        Solid,
        Shell,
        Beam,
        Point,
        Cell,
        MaterialElement
    };
};
}
#endif // ELEMENT_TYPE_F9D2500C_A695_42ec_AE76_94578C5155EA


