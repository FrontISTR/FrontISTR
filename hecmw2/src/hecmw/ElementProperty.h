/* 
 * File:   ElementProperty.h
 * Author: ktakeda
 *
 * Created on 2009/07/16, 12:40
 */

#ifndef _ELEMENTPROPERTY_H_2c62c37d_45fa_4831_adca_c5afc1cd61ae
#define	_ELEMENTPROPERTY_H_2c62c37d_45fa_4831_adca_c5afc1cd61ae

namespace pmw{
struct NumberOfVertex{
    static uint Hexa(){return 8;}    
    static uint Tetra(){return 4;}
    static uint Prism(){return 6;}
    static uint Pyramid(){return 5;}
    static uint Quad(){return 4;}
    static uint Triangle(){return 3;}
    static uint Beam(){return 2;}
    static uint Default(){return 0;}
};
struct NumberOfNode{
    static uint Hexa(){return 8;}
    static uint Hexa2(){return 20;}
    static uint Tetra(){return 4;}
    static uint Tetra2(){return 10;}
    static uint Prism(){return 6;}
    static uint Prism2(){return 15;}
    static uint Pyramid(){return 5;}
    static uint Quad(){return 4;}
    static uint Quad2(){return 8;}
    static uint Triangle(){return 3;}
    static uint Triangle2(){return 6;}
    static uint Beam(){return 2;}
    static uint Beam2(){return 3;}
    static uint Default(){return 0;}
};
struct NumberOfEdge{
    static uint Hexa(){return 12;}
    static uint Tetra(){return 6;}
    static uint Prism(){return 9;}
    static uint Pyramid(){return 8;}
    static uint Quad(){return 4;}
    static uint Triangle(){return 3;}
    static uint Beam(){return 1;}
    static uint Default(){return 0;}
};
struct NumberOfFace{
    static uint Hexa(){return 6;}
    static uint Tetra(){return 4;}
    static uint Prism(){return 5;}
    static uint Pyramid(){return 5;}
    static uint Quad(){return 1;}
    static uint Triangle(){return 1;}
    static uint Beam(){return 0;}
    static uint Default(){return 0;}
};
}

#endif	/* _ELEMENTPROPERTY_H */

