/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ElementProperty.h
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
#ifndef _ELEMENTPROPERTY_H_2c62c37d_45fa_4831_adca_c5afc1cd61ae
#define	_ELEMENTPROPERTY_H_2c62c37d_45fa_4831_adca_c5afc1cd61ae
namespace pmw
{
struct NumberOfVertex {
    static uiint Hexa() {
        return 8;
    }
    static uiint Tetra() {
        return 4;
    }
    static uiint Prism() {
        return 6;
    }
    static uiint Pyramid() {
        return 5;
    }
    static uiint Quad() {
        return 4;
    }
    static uiint Triangle() {
        return 3;
    }
    static uiint Beam() {
        return 2;
    }
    static uiint Point() {
        return 1;
    }
    static uiint Default() {
        return 0;
    }
};
struct NumberOfNode {
    static uiint Hexa() {
        return 8;
    }
    static uiint Hexa2() {
        return 20;
    }
    static uiint Tetra() {
        return 4;
    }
    static uiint Tetra2() {
        return 10;
    }
    static uiint Prism() {
        return 6;
    }
    static uiint Prism2() {
        return 15;
    }
    static uiint Pyramid() {
        return 5;
    }
    static uiint Quad() {
        return 4;
    }
    static uiint Quad2() {
        return 8;
    }
    static uiint Triangle() {
        return 3;
    }
    static uiint Triangle2() {
        return 6;
    }
    static uiint Beam() {
        return 2;
    }
    static uiint Beam2() {
        return 3;
    }
    static uiint Point() {
        return 1;
    }
    static uiint Default() {
        return 0;
    }
};
struct NumberOfEdge {
    static uiint Hexa() {
        return 12;
    }
    static uiint Tetra() {
        return 6;
    }
    static uiint Prism() {
        return 9;
    }
    static uiint Pyramid() {
        return 8;
    }
    static uiint Quad() {
        return 4;
    }
    static uiint Triangle() {
        return 3;
    }
    static uiint Beam() {
        return 1;
    }
    static uiint Point() {
        return 0;
    }
    static uiint Default() {
        return 0;
    }
};
struct NumberOfFace {
    static uiint Hexa() {
        return 6;
    }
    static uiint Tetra() {
        return 4;
    }
    static uiint Prism() {
        return 5;
    }
    static uiint Pyramid() {
        return 5;
    }
    static uiint Quad() {
        return 1;
    }
    static uiint Triangle() {
        return 1;
    }
    static uiint Beam() {
        return 0;
    }
    static uiint Point() {
        return 0;
    }
    static uiint Default() {
        return 0;
    }
};
}
#endif	/* _ELEMENTPROPERTY_H */
