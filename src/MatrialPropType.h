/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   MatrialPropType.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _MATRIAL_Prop_TYPE_H_317b7e89
#define	_MATRIAL_Prop_TYPE_H_317b7e89
namespace pmw{
struct MaterialPropType{
    enum{
        Density,     
        Poisson,     
        YoungModule, 
        Temp_Depend_YoungModule,
        Temp_Depend_Poisson,    
        Linear_Expansion, 
        Thermal_Conductivity, 
        Heat_Transfer_Rate,   
        Cp,
        Cv 
    };
};
}
#endif	/* _MATRIALTYPE_H */
