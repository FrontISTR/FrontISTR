/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ElementProperty.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
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
};
struct NumberOfEdge{
    static uint Hexa(){return 12;}
    static uint Tetra(){return 6;}
    static uint Prism(){return 9;}
    static uint Pyramid(){return 8;}
    static uint Quad(){return 4;}
    static uint Triangle(){return 3;}
    static uint Beam(){return 1;}
};
struct NumberOfFace{
    static uint Hexa(){return 6;}
    static uint Tetra(){return 4;}
    static uint Prism(){return 5;}
    static uint Pyramid(){return 5;}
    static uint Quad(){return 1;}
    static uint Triangle(){return 1;}
    static uint Beam(){return 0;}
};
}
#endif	/* _ELEMENTPROPERTY_H */
