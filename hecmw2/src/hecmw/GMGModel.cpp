/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/GMGModel.cpp
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
#include "GMGModel.h"
using namespace pmw;
CGMGModel::CGMGModel(void)
{
    mpLogger = Utility::CLogger::Instance();
}
CGMGModel::~CGMGModel(void)
{
    for_each(mvAssyModel.begin(), mvAssyModel.end(), DeleteObject());
}
void CGMGModel::initAssyModel(const uiint& num_of_mglevel)
{
    mvAssyModel.resize(num_of_mglevel);
    mNumOfLevel = num_of_mglevel;
}
void CGMGModel::reserveAssyModel(const uiint& num_of_mglevel)
{
    if(num_of_mglevel > mvAssyModel.max_size()) mvAssyModel.reserve(num_of_mglevel);
}
void CGMGModel::setModel(CAssyModel* pAssyModel, const uiint& i)
{
    mvAssyModel[i] = pAssyModel;
}
void CGMGModel::addModel(CAssyModel* pAssyModel, const uiint& i)
{
    if( i < mvAssyModel.size()) {
        mvAssyModel[i]=pAssyModel;
    } else {
        mvAssyModel.push_back(pAssyModel);
    }
}
void CGMGModel::initAssyModel()
{
    mvAssyModel.resize(1);
    mvAssyModel[0]= new CAssyModel;
    mvAssyModel[0]->setMGLevel(0);
    mNumOfLevel = 1;
}
void CGMGModel::resizeAssyModel(const uiint &nNumOfMGLevel)
{
    mvAssyModel.resize(nNumOfMGLevel);
    mNumOfLevel = nNumOfMGLevel;
}
