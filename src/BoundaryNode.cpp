/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BoundaryNode.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "BoundaryNode.h"
using namespace pmw;
CBoundaryNode::CBoundaryNode()
{
    ;
}
CBoundaryNode::~CBoundaryNode()
{
}
void CBoundaryNode::initialize(const uint& bndType,const uint& dof)
{
    mnType = bndType;
    switch(mnType){
        case(BoundaryTypeNode::Disp):
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Accel):
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Load):
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Temp):
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Velo):
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Thermal_Flux):
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        default:
            break;
    }
}
