/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   GMGModel.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
    cout << "~CGMGModel" << endl;
}
void CGMGModel::initAssyModel(const uint& num_of_mglevel)
{ 
    mvAssyModel.resize(num_of_mglevel);
    mNumOfLevel = num_of_mglevel;
    uint num = mvAssyModel.size();
    mpLogger->Info(Utility::LoggerMode::MWDebug, "mvAssyModel.size() == ", num);
}
void CGMGModel::reserveAssyModel(const uint& num_of_mglevel)
{
    if(num_of_mglevel > mvAssyModel.max_size()) mvAssyModel.reserve(num_of_mglevel);
}
void CGMGModel::setModel(CAssyModel* pAssyModel, const uint& i)
{
    mvAssyModel[i] = pAssyModel;
}
void CGMGModel::addModel(CAssyModel* pAssyModel, const uint& i)
{
    if( i < mvAssyModel.size()){
        mvAssyModel[i]=pAssyModel;
    }else{
        mvAssyModel.push_back(pAssyModel);
    }
}
void CGMGModel::resizeAssyModel(const uint &num_of_part)
{
    mvAssyModel.resize(num_of_part);
}
