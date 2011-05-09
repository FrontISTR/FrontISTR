/*
 * File:   Material.h
 * Author: ktakeda
 *
 * 材質データ
 *
 * Created on 2009/05/19, 14:48
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
    uint  mMeshID;//Materialを使用するMeshID => Fileに記述された材質をMeshに割り当てるため.

    map<uint,double,less<uint> > mmValue;//uint==材質PropType => MaterialPropType.h


    //    // ↓ とりあえずイメージとして変数を用意してあるだけなので,削除予定.
    //    double mDensity;
    //    double mPoisson;
    //    double mYoungModule;
    //    double mTempDepend_YoungModule;
    //    double mTempDepend_Poisson;
    //    double mLinear_Expansion;
    //
    //    double mThermalConductivity;//熱伝導
    //    double mCp;//定圧比熱
    //    double mCv;//定積比熱

public:
    // ID
    // --
    void  setID(const uint& id){ mnID = id;}
    uint& getID(){ return mnID;}

    // MeshID
    // --
    void  setMeshID(const uint& mesh_id){ mMeshID= mesh_id;}
    uint& getMeshID(){ return mMeshID;}

    // 名前
    // --
    void setName(std::string& name){ msName= name;}
    std::string& getName(){ return msName;}

    // 材質データ
    // --
    void setValue(const uint& prop_Type, const double& value);
    double& getValue(const uint& prop_Type){ return mmValue[prop_Type];}


    //    // ↓ 削除予定.(具体的なデータは,削除する)
    //    // Material parameter Accessor
    //    //
    //    void setDensity(const double& density){ mDensity=density;}
    //    double& getDensity(){ return mDensity;}
    //
    //    void setPoisson(const double& poisson){ mPoisson=poisson;}
    //    double& getPoisson(){ return mPoisson;}
    //
    //    void setYoungModule(const double& young_module){ mYoungModule=young_module;}
    //    double& getYoungModule(){ return mYoungModule;}
    //
    //    void setTDYoungModule(const double& young_module){ mTempDepend_YoungModule=young_module;}
    //    double& getTDYoungModule(){ return mTempDepend_YoungModule;}
    //
    //    void setTDPoisson(const double& poisson){ mTempDepend_Poisson=poisson;}
    //    double& getTDPoisson(){ return mTempDepend_Poisson;}
    //
    //    void setLinearExpansion(const double& expansion){ mLinear_Expansion=expansion;}
    //    double& getLinearExpansion(){ return mLinear_Expansion;}
    //
    //    void setThermalConductivity(const double& therm_conduct){ mThermalConductivity=therm_conduct;}
    //    double& getThermalConductivity(){ return mThermalConductivity;}
    //
    //    void setCp(const double& val){ mCp=val;}
    //    double& getCp(){return mCp;}
    //
    //    void setCv(const double& val){ mCv=val;}
    //    double& getCv(){return mCv;}

};
}
#endif	/* _MATERIAL_H */

