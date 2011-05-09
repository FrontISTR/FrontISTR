/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   GMGModel.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef GMG_Model_hh_C77868BF_2116_41b3_B667_59BA392FCC5C
#define GMG_Model_hh_C77868BF_2116_41b3_B667_59BA392FCC5C
#include "CommonStd.h"
#include "AssyModel.h"
#include "Material.h"
#include "Logger.h"
namespace pmw{
class CGMGModel
{
public:
static CGMGModel* Instance(){
	static CGMGModel gmgModel;
	return &gmgModel;
}
private:
	CGMGModel(void);
public:
    virtual ~CGMGModel(void);
protected:
    std::vector<CAssyModel*> mvAssyModel;
    uint  mNumOfLevel;
    vector<CMaterial*> mvMaterial;
    Utility::CLogger *mpLogger;
public:
    CAssyModel* getAssyModel(const uint& mgLevel){ return mvAssyModel[mgLevel];}
    uint& getNumOfLevel(){return mNumOfLevel;}
    void  addModel(CAssyModel *pAssyModel, const uint& i);
    void  initAssyModel(const uint& num_of_mglevel);
    void  reserveAssyModel(const uint& num_of_mglevel);
    void  setModel(CAssyModel *pAssyModel, const uint& i);
    void resizeAssyModel(const uint& num_of_part);
    void reserveMaterial(const uint& res_size){ mvMaterial.reserve(res_size);}
    void setMaterial(CMaterial* pMaterial){ mvMaterial.push_back(pMaterial);}
    CMaterial* getMaterial(const uint& i){ return mvMaterial[i];}
    uint getNumOfMaterial(){ return mvMaterial.size();}
};
}
#endif
