/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Material.h
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
#ifndef _MATERIAL_H_878129a5_f446_4878_bf38_8cc90d3fd8a2
#define	_MATERIAL_H_878129a5_f446_4878_bf38_8cc90d3fd8a2
#include <map>
#include "CommonStd.h"
#include "TypeDef.h"
#include "MatrialPropType.h"
namespace pmw
{
class CMaterial
{
public:
    CMaterial();
    virtual ~CMaterial();
private:
    uiint  mnID;
    std::string msName;
    uiint  mMeshID;
    map<uiint,double,less<uiint> > mmValue;
public:
    void  setID(const uiint& id) {
        mnID = id;
    }
    uiint& getID() {
        return mnID;
    }
    void  setMeshID(const uiint& mesh_id) {
        mMeshID= mesh_id;
    }
    uiint& getMeshID() {
        return mMeshID;
    }
    void setName(std::string& name) {
        msName= name;
    }
    std::string& getName() {
        return msName;
    }
    void setValue(const uiint& prop_Type, const double& value);
    double& getValue(const uiint& prop_Type) {
        return mmValue[prop_Type];
    }
};
}
#endif	/* _MATERIAL_H */
