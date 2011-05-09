//
//
//  MW_Str.cpp
//  Structure Analysis "main" Routine
//
//  with pMW static lib
//
//          2009.7.29
//          2009.4.02
//          k.Takeda
#include "./HEC_MW3.h"

using namespace pmw;


int main(int argc, char** argv)
{
    CMW *pMW = CMW::Instance();

    const char* cpath= "./";//pass for Mesh file

    pMW->Initialize(argc, argv, cpath);
    
    pMW->FileRead();

    pMW->Refine();
    pMW->FinalizeRefine();// メモリーの解放

    pMW->FileWrite();

    // FEM start
    uint iAssembleMax, iElemMax, iMeshMax;

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

    // number of linear_algebra_eq, each DOF
    uint nNumOfAlgebra(1);
    uint vDOF[nNumOfAlgebra];
    for(uint iAlge=0; iAlge < nNumOfAlgebra; iAlge++) vDOF[iAlge]= 3;

    //  generate linear_algebra_equation
    pMW->GeneLinearAlgebra(nNumOfAlgebra, vDOF);
    
    //  start of MG-Level loop
    iAssembleMax = pMW->GetNumOfAssembleModel();
    for(uint iAssemble = 0; iAssemble < iAssembleMax; iAssemble++) {
      printf("- Mesh-Level loop ( %d / %d ) \n", iAssemble, iAssembleMax);
      
      // select MG-Level && Algebra-Number
      pMW->SelectAssembleModel( iAssemble );// current MG-Level ( MG-Level==iAssemble )
      pMW->SelectAlgebra(0);                // algebra_eqation num ( number==0 )
      

      //  start of mesh-part loop
      iMeshMax = pMW->GetNumOfMeshPart();
      for(uint iMesh = 0; iMesh < iMeshMax; iMesh++){
	printf("- Mesh loop ( %d / %d ) \n", iMesh, iMeshMax);
	pMW->SelectMeshPart_ID( iMesh );
	
	//  start of element loop
	iElemMax = pMW->getElementSize( iMesh );
	for(uint iElem =0; iElem < iElemMax; iElem++){

          uint numOfLocalNode;
          //double B[6][24], DB[6][24], BDB[24][24], ElemMatrix[24*24], weight, detJ;
          double **B, **DB, **BDB, *ElemMatrix, weight, detJ;

          uint shapeType, elemType;
          pMW->SelectElement_IX(iElem);
          
          elemType = pMW->GetElementType();
          
          uint mat_size;
          // Hexa Element
          if(elemType == pMW->elemtype_hexa() ){
              shapeType = pMW->shapetype_hexa82();
              numOfLocalNode = 8;
              mat_size = numOfLocalNode * 3;// 8Node * 3dof == 24
              B = new double*[6]; DB = new double*[6]; BDB = new double*[mat_size]; ElemMatrix = new double[mat_size*mat_size];
          }
          if(elemType == pMW->elemtype_hexa2()){
              shapeType = pMW->shapetype_hexa202();
              numOfLocalNode = 20;
              mat_size = numOfLocalNode * 3;// 20Node * 3dof == 60
              B = new double*[6]; DB = new double*[6]; BDB = new double*[mat_size]; ElemMatrix = new double[mat_size*mat_size];
          }
          // Tetra Element
          if(elemType == pMW->elemtype_tetra() ){
              shapeType = pMW->shapetype_tetra41();
              numOfLocalNode = 4;
              mat_size = numOfLocalNode * 3;
              B = new double*[6]; DB = new double*[6]; BDB = new double*[mat_size]; ElemMatrix = new double[mat_size*mat_size];
          }
          if(elemType == pMW->elemtype_tetra2()){
              shapeType = pMW->shapetype_tetra104();
              numOfLocalNode = 10;
              mat_size = numOfLocalNode * 3;
              B = new double*[6]; DB = new double*[6]; BDB = new double*[mat_size]; ElemMatrix = new double[mat_size*mat_size];
          }
          // Prism Element
          if(elemType == pMW->elemtype_prism() ){
              shapeType = pMW->shapetype_prism62();
              numOfLocalNode = 6;
              mat_size = numOfLocalNode * 3;
              B = new double*[6]; DB = new double*[6]; BDB = new double*[mat_size]; ElemMatrix = new double[mat_size*mat_size];
          }
          if(elemType == pMW->elemtype_prism2()){
              shapeType = pMW->shapetype_prism159();
              numOfLocalNode = 15;
              mat_size = numOfLocalNode * 3;
              B = new double*[6]; DB = new double*[6]; BDB = new double*[mat_size]; ElemMatrix = new double[mat_size*mat_size];
          }
          uint i;
          for( i=0; i <  6; i++){ B[i] = new double[mat_size]; DB[i] = new double[mat_size];}
          for( i=0; i < mat_size; i++){ BDB[i] = new double[mat_size];}
          
          uint numOfInteg = pMW->NumOfIntegPoint(shapeType);
          
          vvvdouble dNdx;
	  pMW->dNdx(elemType, numOfInteg, iElem, dNdx);

          // 0 clear
	  for(uint i=0; i < mat_size; i++) for(uint j=0; j < mat_size; j++) ElemMatrix[mat_size*i+j] = 0.0;
          
	  // start of Gauss points loop
	  for(uint igauss=0; igauss < numOfInteg; igauss++){
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
            
	    for(uint i=0; i < 6; i++) {
	    for(uint j=0; j < mat_size; j++) {
		DB[i][j] = 0.0;
		for(uint k=0; k < 6; k++)  DB[i][j] += D[i][k]*B[k][j];
	    }
	    }
            
	    for(uint i=0; i < mat_size; i++) {
	    for(uint j=0; j < mat_size; j++) {
		BDB[i][j] = 0.0;
		for(uint k=0; k<6; k++) BDB[i][j] += B[k][i]*DB[k][j];
	    }
	    }
	    pMW->detJacobian(elemType, numOfInteg, igauss, detJ);
	    pMW->Weight(elemType, numOfInteg, 0, weight);
            
	    for(uint i=0; i < mat_size; i++)
            for(uint j=0; j < mat_size; j++)
              ElemMatrix[mat_size*i+j] += BDB[i][j] * weight * detJ;

	  }; // end of Gauss points loop
          
	  pMW->Matrix_Add_Elem(iMesh, iElem, ElemMatrix);
          
	};// end of element loop

        // ------------
        // boundary set
        // ------------
        // 1.boundary face
        uint nNumOfBFMesh = pMW->GetNumOfBoundaryFaceMesh();
        uint iBFMesh;
        for(iBFMesh=0; iBFMesh < nNumOfBFMesh; iBFMesh++) {

            uint nBndType = pMW->GetBNDType_BFaceMesh(iBFMesh);
            uint nBDOF = pMW->GetNumOfDOF_BFaceMesh(iBFMesh);
            uint nNumOfBNode = pMW->GetNumOfBNode_BFaceMesh(iBFMesh);
            uint iBNode;
            for(iBNode=0; iBNode < nNumOfBNode; iBNode++) {
                uint iBNodeID = pMW->GetNodeID_BNode_BFaceMesh(iBFMesh, iBNode);
                // Neumann
                if(nBndType==pMW->getNeumannType()){
                    for(uint idof=0; idof < nBDOF; idof++) { 
                        uint nDOF = pMW->GetDOF_BFaceMesh(iBFMesh,idof);
                        double dBndVal = pMW->GetBNodeValue_BFaceMesh(iBFMesh,iBNode, nDOF, iAssemble);//boundary value
                        
                        pMW->Set_BC_RHS( iMesh, iBNodeID, nDOF, dBndVal);// boundary set for RHS_vector
                    };
                }
                // Dirichlet
                if(nBndType==pMW->getDirichletType()){
                    for(uint idof=0; idof < nBDOF; idof++) {
                        uint nDOF = pMW->GetDOF_BFaceMesh(iBFMesh,idof);
                        double dBndVal = pMW->GetBNodeValue_BFaceMesh(iBFMesh, iBNode, nDOF, iAssemble);//boundary value
                        
                        double dDiagonal = 1.0E+15;        //matrix diagonal value
                        double dRHSVal = dDiagonal*dBndVal;//RHS_vector value
                        pMW->Set_BC_Mat_RHS(iMesh, iBNodeID, nDOF, dDiagonal, dRHSVal);// boundary set for Matrix-Diagonal & RHS_vector
                    };
                }
            };
        };// end of boundary_face_mesh loop

        // 2.boundary edge
        uint nNumOfBEMesh = pMW->GetNumOfBoundaryEdgeMesh();
        uint iBEMesh;
        for(iBEMesh=0; iBEMesh < nNumOfBEMesh; iBEMesh++) {

            uint nBndType = pMW->GetBNDType_BEdgeMesh(iBEMesh);
            uint nBDOF = pMW->GetNumOfDOF_BEdgeMesh(iBEMesh);
            uint nNumOfBNode = pMW->GetNumOfBNode_BEdgeMesh(iBEMesh);
            uint iBNode;
            for(iBNode=0; iBNode < nNumOfBNode; iBNode++) {
                uint iBNodeID = pMW->GetNodeID_BNode_BEdgeMesh(iBEMesh, iBNode);
                // Neumann
                if(nBndType==pMW->getNeumannType()){
                    for(uint idof=0; idof < nBDOF; idof++) {
                        uint nDOF = pMW->GetDOF_BEdgeMesh(iBEMesh,idof);
                        double dBndVal = pMW->GetBNodeValue_BEdgeMesh(iBEMesh,iBNode, nDOF, iAssemble);//boundary value
                        
                        pMW->Set_BC_RHS( iMesh, iBNodeID, nDOF, dBndVal);// boundary set for RHS_vector
                    };
                }
                // Dirichlet
                if(nBndType==pMW->getDirichletType()){
                    for(uint idof=0; idof < nBDOF; idof++) {
                        uint nDOF = pMW->GetDOF_BEdgeMesh(iBEMesh,idof);
                        double dBndVal = pMW->GetBNodeValue_BEdgeMesh(iBEMesh, iBNode, nDOF, iAssemble);//boundary value
                        
                        double dDiagonal = 1.0E+15;        //matrix diagonal value
                        double dRHSVal = dDiagonal*dBndVal;//RHS_vector value
                        pMW->Set_BC_Mat_RHS(iMesh, iBNodeID, nDOF, dDiagonal, dRHSVal);// boundary set for Matrix-Diagonal & RHS_vector
                    };
                }
            };
        };// end of boundary_edge_mesh loop

        // 3.boundary volume
        uint nNumOfBVMesh = pMW->GetNumOfBoundaryVolumeMesh();
        uint iBVMesh;
        for(iBVMesh=0; iBVMesh < nNumOfBVMesh; iBVMesh++) {

            uint nBndType = pMW->GetBNDType_BVolumeMesh(iBVMesh);
            uint nBDOF = pMW->GetNumOfDOF_BVolumeMesh(iBVMesh);
            uint nNumOfBNode = pMW->GetNumOfBNode_BVolumeMesh(iBVMesh);
            uint iBNode;
            for(iBNode=0; iBNode < nNumOfBNode; iBNode++) {
                uint iBNodeID = pMW->GetNodeID_BNode_BVolumeMesh(iBVMesh, iBNode);
                // Neumann
                if(nBndType==pMW->getNeumannType()){
                    for(uint idof=0; idof < nBDOF; idof++) {
                        uint nDOF = pMW->GetDOF_BVolumeMesh(iBVMesh,idof);
                        double dBndVal = pMW->GetBNodeValue_BVolumeMesh(iBVMesh,iBNode, nDOF, iAssemble);//boundary value
                        
                        pMW->Set_BC_RHS( iMesh, iBNodeID, nDOF, dBndVal);// boundary set for RHS_vector
                    };
                }
                // Dirichlet
                if(nBndType==pMW->getDirichletType()){
                    for(uint idof=0; idof < nBDOF; idof++) {
                        uint nDOF = pMW->GetDOF_BVolumeMesh(iBVMesh,idof);
                        double dBndVal = pMW->GetBNodeValue_BVolumeMesh(iBVMesh, iBNode, nDOF, iAssemble);//boundary value
                        
                        double dDiagonal = 1.0E+15;        //matrix diagonal value
                        double dRHSVal = dDiagonal*dBndVal;//RHS_vector value
                        pMW->Set_BC_Mat_RHS(iMesh, iBNodeID, nDOF, dDiagonal, dRHSVal);// boundary set for Matrix-Diagonal & RHS_vector
                    };
                }
            };
        };// end of boundary_edge_mesh loop

        // 4.boundary node
        uint nNumOfBNMesh = pMW->GetNumOfBoundaryNodeMesh();
        uint iBNMesh;
        for(iBNMesh=0; iBNMesh < nNumOfBNMesh; iBNMesh++) {

            uint nBndType = pMW->GetBNDType_BNodeMesh(iBNMesh);
            uint nNumOfBNode = pMW->GetNumOfBNode_BNodeMesh(iBNMesh);
            uint iBNode;
            for(iBNode=0; iBNode < nNumOfBNode; iBNode++) {
                uint nBDOF = pMW->GetNumOfDOF_BNodeMesh(iBNMesh, iBNode);
                uint iBNodeID = pMW->GetNodeID_BNode_BNodeMesh(iBNMesh, iBNode);
                // Neumann
                if(nBndType==pMW->getNeumannType()){
                    for(uint idof=0; idof < nBDOF; idof++) {
                        uint nDOF = pMW->GetDOF_BNodeMesh(iBNMesh, iBNode, idof);
                        double dBndVal = pMW->GetBNodeValue_BNodeMesh(iBNMesh,iBNode,nDOF);//boundary value
                        
                        pMW->Set_BC_RHS( iMesh, iBNodeID, nDOF, dBndVal);// boundary set for RHS_vector
                    };
                }
                // Dirichlet
                if(nBndType==pMW->getDirichletType()){
                    for(uint idof=0; idof < nBDOF; idof++) {
                        uint nDOF = pMW->GetDOF_BNodeMesh(iBNMesh, iBNode, idof);
                        double dBndVal = pMW->GetBNodeValue_BNodeMesh(iBNMesh, iBNode, nDOF);//boundary value
                        
                        double dDiagonal = 1.0E+15;        //matrix diagonal value
                        double dRHSVal = dDiagonal*dBndVal;//RHS_vector value
                        pMW->Set_BC_Mat_RHS(iMesh, iBNodeID, nDOF, dDiagonal, dRHSVal);// boundary set for Matrix-Diagonal & RHS_vector
                    };
                }
            };
        };// end of boundary_node_mesh loop
        
      };// end of mesh-part loop

    };// end of MG-Level loop

    // select MG-Level && Algebra-Number
    pMW->SelectAssembleModel( iAssembleMax - 1 );// finest MG-Level (MG-Level==iAssembleMax-1)
    pMW->SelectAlgebra(0);                       // algebra_eqation num ( number==0 )

    // solve the equation
    uint max_iteration = 100;
    double tolerance = 1.0e-15;
    uint method = 1;       // 1:CG, 2:BiCBSTAB, 3:GPBiCG, 4:GMRES
    uint precondition = 1; // 1:Jacobi, 2:MG, 3:SSOR, 4:ILU
    pMW->Solve( max_iteration, tolerance, method, precondition);
    
    pMW->Finalize();
    return 0;
}

