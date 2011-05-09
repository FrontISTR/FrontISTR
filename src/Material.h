/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Material.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _MATERIAL_H_878129a5_f446_4878_bf38_8cc90d3fd8a2
#define	_MATERIAL_H_878129a5_f446_4878_bf38_8cc90d3fd8a2
#include <map>
#include "CommonStd.h"
#include "TypeDef.h"
#include "MatrialPropType.h"
namespace pmw{
class CMaterial{
public:
    CMaterial();
    virtual ~CMaterial();
private:
    uint  mnID;
    std::string msName;
    uint  mMeshID;
    map<uint,double,less<uint> > mmValue;
public:
    void  setID(const uint& id){ mnID = id;}
    uint& getID(){ return mnID;}
    void  setMeshID(const uint& mesh_id){ mMeshID= mesh_id;}
    uint& getMeshID(){ return mMeshID;}
    void setName(std::string& name){ msName= name;}
    std::string& getName(){ return msName;}
    void setValue(const uint& prop_Type, const double& value);
    double& getValue(const uint& prop_Type){ return mmValue[prop_Type];}
};
}
#endif	/* _MATERIAL_H */
