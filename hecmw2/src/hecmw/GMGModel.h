/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/GMGModel.h
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
#ifndef GMG_Model_hh_C77868BF_2116_41b3_B667_59BA392FCC5C
#define GMG_Model_hh_C77868BF_2116_41b3_B667_59BA392FCC5C
#include "CommonStd.h"
#include "AssyModel.h"
#include "Material.h"
#include "Logger.h"
namespace pmw
{
class CGMGModel
{
public:
    static CGMGModel* Instance() {
        static CGMGModel gmgModel;
        return &gmgModel;
    }
private:
    CGMGModel(void);
public:
    virtual ~CGMGModel(void);
protected:
    std::vector<CAssyModel*> mvAssyModel;
    uiint  mNumOfLevel;
    vector<CMaterial*> mvMaterial;
    Utility::CLogger *mpLogger;
public:
    CAssyModel* getAssyModel(const uiint& mgLevel) {
        return mvAssyModel[mgLevel];
    }
    uiint& getNumOfLevel() {
        return mNumOfLevel;
    }
    void  addModel(CAssyModel *pAssyModel, const uiint& i);
    void  initAssyModel(const uiint& nNumOfMGLevel);
    void  reserveAssyModel(const uiint& nNumOfMGLevel);
    void  setModel(CAssyModel *pAssyModel, const uiint& i);
    void initAssyModel();
    void resizeAssyModel(const uiint& nNumOfMGLevel);
    void reserveMaterial(const uiint& res_size) {
        mvMaterial.reserve(res_size);
    }
    void setMaterial(CMaterial* pMaterial) {
        mvMaterial.push_back(pMaterial);
    }
    CMaterial* getMaterial(const uiint& i) {
        return mvMaterial[i];
    }
    uiint getNumOfMaterial() {
        return mvMaterial.size();
    }
};
}
#endif
