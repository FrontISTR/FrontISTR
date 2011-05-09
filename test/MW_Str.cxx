/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   MW_Str.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/

#include "./HEC_MW3.h"

using namespace pmw;


int main(int argc, char** argv)
{
    CMW *pMW = CMW::Instance();

    const char* cpath= "./";//pass of mw3.cnt

    pMW->Initialize(argc, argv, cpath);
    
    pMW->FileRead();

    pMW->Refine();

    pMW->FileWrite();

    // FEM start
    uint iAssembleMax, iElemMax, iMeshMax, numOfInteg;
    const uint ishape = ShapeType::Hexa82;
    const uint numOfLocalNode = 8;
    vvvdouble dNdx;
    double B[6][24], DB[6][24], BDB[24][24], ElemMatrix[24*24], weight, J;
	
    // Material Properties
    double nyu = 0.3;
    double E = 1.0e9;
    double coef = ((1.0-nyu)*E)/((1.0+nyu)*(1.0-2*nyu));
    double D[6][6] = { {coef, coef*nyu/(1.0-nyu), coef*nyu/(1.0-nyu), 0, 0, 0},
		       {coef*nyu/(1.0-nyu), coef, coef*nyu/(1.0-nyu), 0, 0, 0},
		       {coef*nyu/(1.0-nyu), coef*nyu/(1.0-nyu), coef, 0, 0, 0},
		       {0, 0, 0, coef*(1-2*nyu)/(2*(1.0-nyu)), 0, 0},
		       {0, 0, 0, 0, coef*(1-2*nyu)/(2*(1.0-nyu)), 0},
		       {0, 0, 0, 0, 0, coef*(1-2*nyu)/(2*(1.0-nyu))} };

    //  start of mesh-level loop
    pMW->GetNumOfAssembleModel(iAssembleMax);
    for(uint iAssemble = 0; iAssemble < iAssembleMax; iAssemble++) {
      printf("- Mesh-Level loop ( %d / %d ) \n", iAssemble, iAssembleMax);
      pMW->SelectAssembleModel( iAssemble );
      pMW->Initialize_Matrix();
      pMW->Initialize_Vector();

      //  start of mesh-part loop
      pMW->GetNumOfMeshPart(iMeshMax);
      for(uint iMesh = 0; iMesh < iMeshMax; iMesh++){
	printf("- Mesh loop ( %d / %d ) \n", iMesh, iMeshMax);
	pMW->SelectMeshPart_ID( iMesh );
	pMW->NumOfIntegPoint(ishape, numOfInteg);

	//  start of element loop
	iElemMax = pMW->getElementSize( iMesh );
	for(uint iElem =0; iElem < iElemMax; iElem++){
	  pMW->dNdx(ishape, numOfInteg, iElem, dNdx);
	  for(uint i=0; i<24; i++) for(uint j=0; j<24; j++)
				  ElemMatrix[24*i+j] = 0.0;
	  // start of Gauss points loop
	  for(uint igauss=0; igauss<numOfInteg; igauss++){
	    for(uint iLocalNode=0; iLocalNode < numOfLocalNode; iLocalNode++){
	      B[0][iLocalNode*3   ] = dNdx[igauss][iLocalNode][0];//dNdX
	      B[0][iLocalNode*3 +1] = 0.0;
	      B[0][iLocalNode*3 +2] = 0.0;
	      B[1][iLocalNode*3   ] = 0.0;
	      B[1][iLocalNode*3 +1] = dNdx[igauss][iLocalNode][1];//dNdY
	      B[1][iLocalNode*3 +2] = 0.0;
	      B[2][iLocalNode*3   ] = 0.0;
	      B[2][iLocalNode*3 +1] = 0.0;
	      B[2][iLocalNode*3 +2] = dNdx[igauss][iLocalNode][2];//dNdZ
	      B[3][iLocalNode*3   ] = dNdx[igauss][iLocalNode][1];//dNdY
	      B[3][iLocalNode*3 +1] = dNdx[igauss][iLocalNode][0];//dNdX
	      B[3][iLocalNode*3 +2] = 0.0;
	      B[4][iLocalNode*3   ] = 0.0;
	      B[4][iLocalNode*3 +1] = dNdx[igauss][iLocalNode][2];//dNdZ
	      B[4][iLocalNode*3 +2] = dNdx[igauss][iLocalNode][1];//dNdY
	      B[5][iLocalNode*3   ] = dNdx[igauss][iLocalNode][2];//dNdZ
	      B[5][iLocalNode*3 +1] = 0.0;
	      B[5][iLocalNode*3 +2] = dNdx[igauss][iLocalNode][0];//dNdX
	    }
	    for(uint i=0; i<6; i++) {
	      for(uint j=0; j<24; j++) {
		DB[i][j] = 0.0;
		for(uint k=0; k<6; k++) {
		  DB[i][j] += D[i][k]*B[k][j];
		}
	      }
	    }
	    for(uint i=0; i<24; i++) {
	      for(uint j=0; j<24; j++) {
		BDB[i][j] = 0.0;
		for(uint k=0; k<6; k++) {
		  BDB[i][j] += B[k][i]*DB[k][j];
		}
	      }
	    }
	    pMW->detJacobian(ishape, numOfInteg, igauss, J);
	    pMW->Weight(ishape, numOfInteg, 0, weight);
	    for(uint i=0; i<24; i++) for(uint j=0; j<24; j++)
				  ElemMatrix[24*i+j] += BDB[i][j] * weight * J;
	  }; // end of Gauss points loop
	  
	  pMW->Matrix_Add_Elem(iMesh, iElem, ElemMatrix);
	}; // end of element loop
	
	pMW->Sample_Set_BC( iMesh );
      }; // end of mesh loop
      
      pMW->StoreMatrix();
    }; // end of Mesh-Level loop
    
    // load the equation for the finest mesh (iAssemble=iAssembleMax-1)
    pMW->SelectAssembleModel( iAssembleMax - 1 );
    pMW->LoadMatrix();

    // solve the equation
    uint max_iteration = 1000;
    double tolerance = 1.0e-4;
    uint method = 1;       // 1:CG, 2:BiCBSTAB, 3:GPBiCG, 4:GMRES
    uint precondition = 1; // 1:Jacobi, 2:MG, 3:SSOR, 4:ILU
    pMW->Solve( max_iteration, tolerance, method, precondition);
    
    pMW->Finalize();
    return 0;
}

