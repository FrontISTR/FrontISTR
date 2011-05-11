//
//	GMGModel.h
//
//				2009.05.01
//				2008.11.10
//				k.Takeda
#ifndef GMG_Model_hh_C77868BF_2116_41b3_B667_59BA392FCC5C
#define GMG_Model_hh_C77868BF_2116_41b3_B667_59BA392FCC5C

#include "CommonStd.h"
#include "AssyModel.h"

#include "Material.h"

#include "Logger.h"

namespace pmw{
class CGMGModel
{
////public:
////	static CGMGModel* Instance(){
////		if(!pInstance)
////		    pInstance = new CGMGModel();
////		return pInstance;
////	}
////private:
////	static CGMGModel* pInstance;
////	CGMGModel(void);
public:
static CGMGModel* Instance(){
	static CGMGModel gmgModel;
	return &gmgModel;
}
private:
	CGMGModel(void);
public:
    //CGMGModel(void);
    virtual ~CGMGModel(void);


protected:
    std::vector<CAssyModel*> mvAssyModel;
    uiint  mNumOfLevel;//MultiGrid Level数(階層数)

    vector<CMaterial*> mvMaterial;//材質数だけのMaterialオブジェクト

    Utility::CLogger *mpLogger;

public:
    // AssemblyModel ( MultiGrid Level )
    CAssyModel* getAssyModel(const uiint& mgLevel){ return mvAssyModel[mgLevel];}

    // MultiGrid Level
    //void  setNumOfLevel(const uint& num_of_level){ mNumOfLevel= num_of_level;}
    uiint& getNumOfLevel(){return mNumOfLevel;}


    // GMG
    void  addModel(CAssyModel *pAssyModel, const uiint& i);
    void  initAssyModel(const uiint& nNumOfMGLevel);      // FileReaderRefine -> Factory::GeneAssyModelから使用  
    void  reserveAssyModel(const uiint& nNumOfMGLevel);   // FileReaderRefine -> Factory::GeneAssyModelから使用
    void  setModel(CAssyModel *pAssyModel, const uiint& i);
    
    void initAssyModel();                                // CMW から初期化する場合の関数
    void resizeAssyModel(const uiint& nNumOfMGLevel);    // CMW::Refine -> Factory::GeneAssyModelから使用



    // Material
    void reserveMaterial(const uiint& res_size){ mvMaterial.reserve(res_size);}
    void setMaterial(CMaterial* pMaterial){ mvMaterial.push_back(pMaterial);}
    CMaterial* getMaterial(const uiint& i){ return mvMaterial[i];}
    uiint getNumOfMaterial(){ return mvMaterial.size();}

};
}
#endif
