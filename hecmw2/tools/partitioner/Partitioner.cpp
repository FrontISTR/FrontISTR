/*
 ----------------------------------------------------------
|
| Software Name :part Ver 0.1 beta
|
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Partitioner.h"
#include <limits.h>

CPartitioner::CPartitioner(CMW* mw, uiint n)
{
    iMeshMax = 1;
    nparts = n;
    pMW = mw;
    nn = 0;
    ne = 0;
    numGlobalCommMesh = 0;
}
CPartitioner::~CPartitioner()
{
    free(elmnts);
    free(epart);
    free(npart);
}

void CPartitioner::prepare()
{
    // prepare
    // count number of nodes and elements
    //uiint iMeshMax = pMW->GetNumOfMeshPart();
    for(uiint iMesh = 0; iMesh < iMeshMax; iMesh++) {
        pMW->SelectMeshPart_IX(iMesh);

        uiint iElemMax = pMW->getElementSize(iMesh);
        nn += pMW->getNodeSize(iMesh);
        for(uiint iElem=0; iElem < iElemMax; iElem++) {
            pMW->SelectElement_IX(iElem);

            uiint elemType = pMW->GetElementType();
            if(elemType == pMW->elemtype_tetra()) {
                ne ++;
            } else {
                //TODO
            }
        }
    }

    elmnts = (idxtype*)malloc(sizeof(idxtype) * ne * 4);
    epart = (idxtype*)malloc(sizeof(idxtype) * ne);
    npart = (idxtype*)malloc(sizeof(idxtype) * nn);

    // prepare mesh data for metis
    int count = 0;
    for(uiint iMesh = 0; iMesh < iMeshMax; iMesh++) {
        pMW->SelectMeshPart_IX(iMesh);

        uiint iElemMax = pMW->getElementSize(iMesh);
        for(uiint iElem=0; iElem < iElemMax; iElem++) {
            pMW->SelectElement_IX(iElem);

            uiint elemType = pMW->GetElementType();
            if(elemType == pMW->elemtype_tetra()) {
                iint* nodeIDArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                pMW->GetElementVertNodeIndex(nodeIDArray);
                for(uiint i=0; i<pMW->GetNumOfElementVert(); i++) {
                    elmnts[count] = nodeIDArray[i];
                    count++;
                }
            } else {
                //TODO
            }
        }
    }

    etype = 2;   // tetra only in this version
    numflag = 0; // arrray start index (0 or 1)

}
void CPartitioner::partition()
{

    //partitioning
    METIS_PartMeshNodal(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);

    //post process
    globalComm.assign((nparts * nparts - nparts)/2, 0);
    for(uiint iRegion=0; iRegion < nparts; iRegion++) {
        int index = 0;
        for(uiint iMesh = 0; iMesh < iMeshMax; iMesh++) {
            pMW->SelectMeshPart_IX(iMesh);
            uiint iElemMax = pMW->getElementSize(iMesh);
            for(uiint iElem=0; iElem < iElemMax; iElem++) {
                pMW->SelectElement_IX(iElem);
                uiint elemType = pMW->GetElementType();
                if(epart[index] == iRegion) {
                    //CommFace
                    for(uiint i=0; i<4; i++) { //TODO GetNumOfElementFace should be implemented in MW.
                        if(pMW->GetElementFaceElementID(i)) {
                            uiint elementID = pMW->GetElementFaceElementID(i);
                            uiint elementIndex = pMW->getElementIndex(pMW->GetElementFaceElementID(i));
                            if(epart[elementIndex] > iRegion) {
                                globalComm.at(getCommMeshID(iRegion, epart[elementIndex], nparts))++;
                            }
                        }
                    }
                    //CommEdge;
                    for(uiint i=0; i<6; i++) {
                        for(uiint j=0; j<pMW->GetNumOfElementEdgeElement(i); j++) {
                            uiint elementID = pMW->GetElementEdgeElementID(i,j);
                            uiint elementIndex = pMW->getElementIndex(elementID);
                            if(epart[elementIndex] > iRegion) {
                                //printf("region = %d index = %d, edgeIndex = %d, neighborElementId = %d neighborELementIndex = %d neighborRegion = %d\n", iRegion, index, i, elementID, elementIndex, epart[elementIndex]);
                                globalComm.at(getCommMeshID(iRegion, epart[elementIndex], nparts))++;
                            }
                        }
                    }
                    //CommPoint;
                    iint* nodeIndexArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                    pMW->GetElementVertNodeIndex(nodeIndexArray);
                    for(uiint i=0; i<4; i++) {
                        uiint nodeID = pMW->getNodeID(nodeIndexArray[i]);
                        uiint numOfAggregateElement = pMW->getNumOfAggregateElement(nodeID);
                        for(uiint j=0; j<numOfAggregateElement; j++) {
                            uiint elementID = pMW->getAggregateElementID(nodeID, j);
                            uiint elementIndex = pMW->getElementIndex(elementID);
                            uiint elementRegion = epart[elementIndex];
                            if(elementRegion > iRegion) {
                                globalComm.at(getCommMeshID(iRegion, epart[elementIndex], nparts))++;
                            }
                        }
                    }
                    free(nodeIndexArray);
                }
                index++;
            }
        }
    }
    numGlobalCommMesh = 0;
    for(uiint i=0; i<(nparts * nparts - nparts)/2; i++) {
        if(globalComm.at(i) > 0) {
            numGlobalCommMesh++;
        }
    }

    /*
    FILE* fp = fopen("epart.txt", "w");
    for(int i = 0; i < ne; i++){
        fprintf(fp, "%d %d %d\n", i, pMW->getElementID(i), epart[i]);
    }
    fclose(fp);
    */
}

void CPartitioner::post()
{
    //make nodes and elements table;
    npartTable.assign(nparts, vector<uiint>(nn, -1));
    epartTable.assign(nparts, vector<uiint>(ne, -1));

    for(uiint iRegion=0; iRegion < nparts; iRegion++) {

        int index = 0;
        for(uiint iMesh = 0; iMesh < iMeshMax; iMesh++) {
            pMW->SelectMeshPart_IX(iMesh);
            //elem loop
            uiint iElemMax = pMW->getElementSize(iMesh);
            for(uiint iElem=0; iElem < iElemMax; iElem++) {
                //pMW->LoggerInfo(nInfo, "%s%d%s%d", "-- Element loop", iElem, "MaxElem:", iElemMax-1);
                pMW->SelectElement_IX(iElem);

                uiint elemType = pMW->GetElementType();
                if(elemType == pMW->elemtype_tetra()) {
                    iint* nodeIndexArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                    pMW->GetElementVertNodeIndex(nodeIndexArray);

                    if(epart[index] == iRegion) {
                        epartTable.at(iRegion).at(index) = epart[index];
                        //printf("element index = %d id = %d region = %d\n", index, pMW->getElementID(index), epartTmp.at(index));
                        for(uiint i=0; i<pMW->GetNumOfElementVert(); i++) {
                            npartTable.at(iRegion).at(nodeIndexArray[i]) = iRegion;
                        }
                    }
                }
                index++;
            }
        }
    }

    //count nodes
    numNodesTable.assign(nparts, 0);
    for(uiint iRegion=0; iRegion < nparts; iRegion++) {
        uiint nodeCount = 0;
        for(uiint i=0; i<npartTable.at(iRegion).size(); i++) {
            if(npartTable.at(iRegion).at(i) == iRegion) {
                nodeCount++;
            }
        }
        numNodesTable.at(iRegion) = nodeCount;
    }

    //count elements
    numElementsTable.assign(nparts, 0);
    for(uiint iRegion=0; iRegion < nparts; iRegion++) {
        uiint elementCount = 0;
        for(uiint i=0; i<epartTable.at(iRegion).size(); i++) {
            if(epartTable.at(iRegion).at(i) == iRegion) {
                elementCount++;
            }
        }
        numElementsTable.at(iRegion) = elementCount;
    }

    //initialize
    commFaceNumTable.assign(nparts, vector<uiint>(nparts, 0));
    commNodeTable.assign(nparts, vector<vector<uiint> >(nparts, vector<uiint>(nn, -1)));
    commNodeIndexTable.assign(nparts, vector<vector<uiint> >(nparts, vector<uiint>(nn, -1)));
    commFaceElementTable.assign(nparts, vector<uiint>(0));
    commFaceIndexTable.assign(nparts, vector<uiint>(0));
    commFaceElementVertTable.assign(nparts, vector<uiint>(0));
    commFaceIndexVertTable.assign(nparts, vector<uiint>(0));
    commFaceRegionTable.assign(nparts, vector<uiint>(0));
    commFaceTypeTable.assign(nparts, vector<uiint>(0));

    //make commFace
    for(uiint iRegion=0; iRegion < nparts; iRegion++) {
        //printf("region=%d\n", iRegion);

        int index = 0;
        for(uiint iMesh = 0; iMesh < iMeshMax; iMesh++) {
            //pMW->LoggerInfo(nInfo, "%s%d%s%d", "-- Mesh loop", iMesh, "MaxMeshParts:", iMeshMax-1);
            pMW->SelectMeshPart_IX(iMesh);

            //elem loop
            uiint iElemMax = pMW->getElementSize(iMesh);
            for(uiint iElem=0; iElem < iElemMax; iElem++) {
                //pMW->LoggerInfo(nInfo, "%s%d%s%d", "-- Element loop", iElem, "MaxElem:", iElemMax-1);
                pMW->SelectElement_IX(iElem);

                uiint elemType = pMW->GetElementType();
                if(elemType == pMW->elemtype_tetra()) {
                    iint* nodeIndexArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                    pMW->GetElementVertNodeIndex(nodeIndexArray);

                    if(epart[index] == iRegion) {

                        for(uiint i=0; i<4; i++) {
                            if(pMW->GetElementFaceElementID(i)) {
                                uiint elementID = pMW->GetElementFaceElementID(i);
                                uiint elementIndex = pMW->getElementIndex(pMW->GetElementFaceElementID(i));
                                if(epart[elementIndex] > iRegion) {
                                    //printf("face:%d elementID:%d elementIndex:%d region:%d\n", i, elementID, elementIndex, epart[elementIndex]);
                                    commFaceNumTable.at(iRegion).at(epart[elementIndex])++;
                                    //commFaceElement.push_back(elementIndex);
                                    commFaceElementTable.at(iRegion).push_back(index);
                                    commFaceIndexTable.at(iRegion).push_back(i);
                                    commFaceRegionTable.at(iRegion).push_back(epart[elementIndex]);
                                    commFaceTypeTable.at(iRegion).push_back(3);//triangle

                                    if(i==0) {
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[0]) = epart[elementIndex];
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[1]) = epart[elementIndex];
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[2]) = epart[elementIndex];
                                    } else if(i==1) {
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[0]) = epart[elementIndex];
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[3]) = epart[elementIndex];
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[1]) = epart[elementIndex];
                                    } else if(i==2) {
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[1]) = epart[elementIndex];
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[3]) = epart[elementIndex];
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[2]) = epart[elementIndex];
                                    } else if(i==3) {
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[0]) = epart[elementIndex];
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[2]) = epart[elementIndex];
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndexArray[3]) = epart[elementIndex];
                                    }

                                    //相手要素における面番号の検索
                                    pMW->SelectElement_IX(elementIndex);
                                    uiint faceNum = -1;
                                    for(uiint j=0; j<4; j++) {
                                        if(pMW->GetElementFaceElementID(j) == pMW->getElementID(iElem)) {
                                            faceNum = j;
                                            break;
                                        }
                                    }
                                    pMW->SelectElement_IX(iElem);

                                    //相手要素の登録
                                    commFaceElementVertTable.at(epart[elementIndex]).push_back(index);
                                    commFaceIndexVertTable.at(epart[elementIndex]).push_back(i);
                                    commFaceElementVertTable.at(iRegion).push_back(elementIndex);
                                    commFaceIndexVertTable.at(iRegion).push_back(index);

                                    commFaceNumTable.at(epart[elementIndex]).at(iRegion)++;
                                    commFaceElementTable.at(epart[elementIndex]).push_back(elementIndex);
                                    commFaceIndexTable.at(epart[elementIndex]).push_back(faceNum);
                                    commFaceRegionTable.at(epart[elementIndex]).push_back(iRegion);
                                    commFaceTypeTable.at(epart[elementIndex]).push_back(3);//triangle
                                    if(i==0) {
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[0]) = iRegion;
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[1]) = iRegion;
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[2]) = iRegion;
                                    } else if(i==1) {
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[0]) = iRegion;
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[3]) = iRegion;
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[1]) = iRegion;
                                    } else if(i==2) {
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[1]) = iRegion;
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[3]) = iRegion;
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[2]) = iRegion;
                                    } else if(i==3) {
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[0]) = iRegion;
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[2]) = iRegion;
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndexArray[3]) = iRegion;
                                    }


                                }
                            }
                        }
                    }

                    index++;
                } else {
                }
            }
        }

        index = 0;
        for(uiint iMesh = 0; iMesh < iMeshMax; iMesh++) {
            //pMW->LoggerInfo(nInfo, "%s%d%s%d", "-- Mesh loop", iMesh, "MaxMeshParts:", iMeshMax-1);
            pMW->SelectMeshPart_IX(iMesh);

            //elem loop
            uiint iElemMax = pMW->getElementSize(iMesh);
            for(uiint iElem=0; iElem < iElemMax; iElem++) {
                //pMW->LoggerInfo(nInfo, "%s%d%s%d", "-- Element loop", iElem, "MaxElem:", iElemMax-1);
                pMW->SelectElement_IX(iElem);

                uiint elemType = pMW->GetElementType();
                if(elemType == pMW->elemtype_tetra()) {
                    iint* nodeIndexArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                    pMW->GetElementVertNodeIndex(nodeIndexArray);

                    if(epart[index] == iRegion) {

                        for(uiint i=0; i<6; i++) {
                            for(uiint j=0; j<pMW->GetNumOfElementEdgeElement(i); j++) {
                                uiint elementID = pMW->GetElementEdgeElementID(i,j);
                                uiint elementIndex = pMW->getElementIndex(elementID);

                                uiint nodeIndex1;
                                uiint nodeIndex2;
                                if(i==0) {
                                    nodeIndex1 = nodeIndexArray[0];
                                    nodeIndex2 = nodeIndexArray[1];
                                } else if(i==1) {
                                    nodeIndex1 = nodeIndexArray[1];
                                    nodeIndex2 = nodeIndexArray[2];
                                } else if(i==2) {
                                    nodeIndex1 = nodeIndexArray[2];
                                    nodeIndex2 = nodeIndexArray[0];
                                } else if(i==3) {
                                    nodeIndex1 = nodeIndexArray[0];
                                    nodeIndex2 = nodeIndexArray[3];
                                } else if(i==4) {
                                    nodeIndex1 = nodeIndexArray[1];
                                    nodeIndex2 = nodeIndexArray[3];
                                } else if(i==5) {
                                    nodeIndex1 = nodeIndexArray[2];
                                    nodeIndex2 = nodeIndexArray[3];
                                }

                                //相手要素における辺番号の検索
                                pMW->SelectElement_IX(elementIndex);
                                iint* nodeIndexArrayVs = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                                pMW->GetElementVertNodeIndex(nodeIndexArrayVs);
                                uiint edgeNum;
                                for(uiint j=0; j<6; j++) {
                                    uiint nodeIndexVs1;
                                    uiint nodeIndexVs2;
                                    if(j==0) {
                                        nodeIndexVs1 = nodeIndexArrayVs[0];
                                        nodeIndexVs2 = nodeIndexArrayVs[1];
                                    } else if(j==1) {
                                        nodeIndexVs1 = nodeIndexArrayVs[1];
                                        nodeIndexVs2 = nodeIndexArrayVs[2];
                                    } else if(j==2) {
                                        nodeIndexVs1 = nodeIndexArrayVs[2];
                                        nodeIndexVs2 = nodeIndexArrayVs[0];
                                    } else if(j==3) {
                                        nodeIndexVs1 = nodeIndexArrayVs[0];
                                        nodeIndexVs2 = nodeIndexArrayVs[3];
                                    } else if(j==4) {
                                        nodeIndexVs1 = nodeIndexArrayVs[1];
                                        nodeIndexVs2 = nodeIndexArrayVs[3];
                                    } else if(j==5) {
                                        nodeIndexVs1 = nodeIndexArrayVs[2];
                                        nodeIndexVs2 = nodeIndexArrayVs[3];
                                    }

                                    if(nodeIndex1 == nodeIndexVs1 && nodeIndex2 == nodeIndexVs2) {
                                        edgeNum = j;
                                        break;
                                    }
                                    if(nodeIndex1 == nodeIndexVs2 && nodeIndex2 == nodeIndexVs1) {
                                        edgeNum = j;
                                        break;
                                    }
                                }
                                pMW->SelectElement_IX(iElem);

                                //登録
                                if(epart[elementIndex] > iRegion) {
                                    uiint t1 = commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndex1);
                                    uiint t2 = commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndex2);
                                    //printf("region = %d index = %d, edgeIndex = %d, neighborElementId = %d neighborELementIndex = %d neighborRegion = %d  %d-%d %d %d\n", iRegion, index, i, elementID, elementIndex, epart[elementIndex], nodeIndex1, nodeIndex2, t1, t2);
                                    if(t1 == -1) {
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndex1) = epart[elementIndex];
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndex1) = iRegion;
                                    }
                                    if(t2 == -1) {
                                        commNodeTable.at(iRegion).at(epart[elementIndex]).at(nodeIndex2) = epart[elementIndex];
                                        commNodeTable.at(epart[elementIndex]).at(iRegion).at(nodeIndex2) = iRegion;
                                    }
                                    if(t1 == -1 || t2 == -1) {
                                        //CommFaceにbeamを追加
                                        commFaceNumTable.at(iRegion).at(epart[elementIndex])++;
                                        commFaceElementTable.at(iRegion).push_back(index);
                                        commFaceIndexTable.at(iRegion).push_back(i);//辺番号
                                        commFaceRegionTable.at(iRegion).push_back(epart[elementIndex]);
                                        commFaceTypeTable.at(iRegion).push_back(1);//beam

                                        //相手
                                        commFaceNumTable.at(epart[elementIndex]).at(iRegion)++;
                                        commFaceElementTable.at(epart[elementIndex]).push_back(elementIndex);
                                        commFaceIndexTable.at(epart[elementIndex]).push_back(edgeNum);//辺番号
                                        commFaceRegionTable.at(epart[elementIndex]).push_back(iRegion);
                                        commFaceTypeTable.at(epart[elementIndex]).push_back(1);//beam

                                        commFaceElementVertTable.at(epart[elementIndex]).push_back(index);
                                        commFaceIndexVertTable.at(epart[elementIndex]).push_back(i);
                                        commFaceElementVertTable.at(iRegion).push_back(elementIndex);
                                        commFaceIndexVertTable.at(iRegion).push_back(index);


                                    }
                                } else if(epart[elementIndex] < iRegion) {

                                }
                            }
                        }
                    }

                    index++;
                } else {
                }
            }
        }


        //search point
        index = 0;
        for(uiint iMesh = 0; iMesh < iMeshMax; iMesh++) {
            //pMW->LoggerInfo(nInfo, "%s%d%s%d", "-- Mesh loop", iMesh, "MaxMeshParts:", iMeshMax-1);
            pMW->SelectMeshPart_IX(iMesh);

            //elem loop
            uiint iElemMax = pMW->getElementSize(iMesh);
            for(uiint iElem=0; iElem < iElemMax; iElem++) {
                //pMW->LoggerInfo(nInfo, "%s%d%s%d", "-- Element loop", iElem, "MaxElem:", iElemMax-1);
                pMW->SelectElement_IX(iElem);

                uiint elemType = pMW->GetElementType();
                if(elemType == pMW->elemtype_tetra()) {
                    iint* nodeIndexArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                    pMW->GetElementVertNodeIndex(nodeIndexArray);

                    if(epart[index] == iRegion) {

                        for(uiint i=0; i<4; i++) {
                            uiint nodeID = pMW->getNodeID(nodeIndexArray[i]);
                            uiint numOfAggregateElement = pMW->getNumOfAggregateElement(nodeID);
                            //printf("node ID = %d num = %d\n", nodeID, numOfAggregateElement);
                            for(uiint j=0; j<numOfAggregateElement; j++) {
                                uiint elementID = pMW->getAggregateElementID(nodeID, j);
                                uiint elementIndex = pMW->getElementIndex(elementID);
                                uiint elementRegion = epart[elementIndex];

                                //相手要素における点番号の検索
                                pMW->SelectElement_IX(elementIndex);
                                iint* nodeIndexArrayVs = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                                pMW->GetElementVertNodeIndex(nodeIndexArrayVs);
                                uiint nodeNum;
                                for(uiint j=0; j<pMW->GetNumOfElementVert(); j++) {
                                    if(nodeIndexArrayVs[j] == nodeIndexArray[i]) {
                                        nodeNum = j;
                                        break;
                                    }
                                }
                                pMW->SelectElement_IX(iElem);

                                uiint t1 = commNodeTable.at(iRegion).at(elementRegion).at(nodeIndexArray[i]);
                                if(elementRegion > iRegion && t1 == -1) {
                                    //printf("aggregate element ID = %d index = %d region = %d\n", elementID, elementIndex, elementRegion);
                                    //CommFaceにpointを追加
                                    commNodeTable.at(iRegion).at(elementRegion).at(nodeIndexArray[i]) = elementRegion;;
                                    commFaceNumTable.at(iRegion).at(elementRegion)++;
                                    commFaceElementTable.at(iRegion).push_back(index);
                                    commFaceIndexTable.at(iRegion).push_back(i);//接点番号
                                    commFaceRegionTable.at(iRegion).push_back(elementRegion);
                                    commFaceTypeTable.at(iRegion).push_back(0);//point

                                    commNodeTable.at(elementRegion).at(iRegion).at(nodeIndexArray[i]) = iRegion;;
                                    commFaceNumTable.at(elementRegion).at(iRegion)++;
                                    commFaceElementTable.at(elementRegion).push_back(elementIndex);
                                    commFaceIndexTable.at(elementRegion).push_back(nodeNum);//接点番号
                                    commFaceRegionTable.at(elementRegion).push_back(iRegion);
                                    commFaceTypeTable.at(elementRegion).push_back(0);//point

                                    commFaceElementVertTable.at(epart[elementIndex]).push_back(index);
                                    commFaceIndexVertTable.at(epart[elementIndex]).push_back(i);
                                    commFaceElementVertTable.at(iRegion).push_back(elementIndex);
                                    commFaceIndexVertTable.at(iRegion).push_back(index);

                                }
                            }
                        }
                    }

                    index++;
                } else {
                }
            }
        }

        for(int i=0; i<nparts; i++) {
            if(commFaceNumTable.at(iRegion).at(i) > 0) {
                uiint index = 0;
                for(int j=0; j<commNodeTable.at(iRegion).at(i).size(); j++) {
                    if(commNodeTable.at(iRegion).at(i).at(j) == i) {
                        commNodeIndexTable.at(iRegion).at(i).at(j) = index;
                        index++;
                    }
                }
            }
        }
    }

}

void CPartitioner::write(string prefix)
{
    //post process + write
    for(uiint iRegion=0; iRegion < nparts; iRegion++) {
        //FILE* fp=stdout;
        char filename[255];;
        sprintf(filename, "%s.%d.msh", prefix.c_str(), iRegion);
        cout << "writing " << filename << endl;
        FILE* fp = fopen(filename, "w");

        //write assuming iMeshMax=0 TODO
        fprintf(fp, "AssyModel\n");
        fprintf(fp, "1 0 0\n");
        fprintf(fp, "0 0\n");
        fprintf(fp, "%d\n", numGlobalCommMesh); //num of globall comm
        for(uiint myrank=0; myrank < nparts; myrank++) {
            for(uiint rank=0; rank < nparts; rank++) {
                if(myrank < rank && globalComm.at(getCommMeshID(myrank, rank, nparts)) > 0) {
                    fprintf(fp, "0 %d %d\n", myrank, rank);
                }
            }
        }
        fprintf(fp, "End\n\n");

        fprintf(fp, "Node\n");
        uiint nodeCount = 0;
        uiint min = INT_MAX;
        uiint max = 0;
        for(uiint i=0; i<npartTable.at(iRegion).size(); i++) {
            if(npartTable.at(iRegion).at(i) == iRegion) {
                nodeCount++;
                if(pMW->getNodeID(i) < min) {
                    min = pMW->getNodeID(i);
                }
                if(pMW->getNodeID(i) > max) {
                    max = pMW->getNodeID(i);
                }
            }
        }
        fprintf(fp, "%d 0 %d %d\n", nodeCount, max, min);
        for(uiint i=0; i<npartTable.at(iRegion).size(); i++) {
            if(npartTable.at(iRegion).at(i) == iRegion) {
                uiint id = pMW->getNodeID(i);
                double x,y,z;
                pMW->GetNodeCoord(id, x, y, z);
                fprintf(fp, "V %d %d %d %f %f %f\n", pMW->GetNumOfScalar(id), pMW->GetNumOfVector(id), id, x, y, z);
            }
        }
        fprintf(fp, "End\n\n");


        fprintf(fp, "Element\n");
        uiint elementCount = 0;
        min = INT_MAX;
        max = 0;
        for(uiint i=0; i<epartTable.at(iRegion).size(); i++) {
            if(epartTable.at(iRegion).at(i) == iRegion) {
                elementCount++;
                if(pMW->getElementID(i) < min) {
                    min = pMW->getElementID(i);
                }
                if(pMW->getElementID(i) > max) {
                    max = pMW->getElementID(i);
                }
            }
        }
        fprintf(fp, "%d 0 %d %d\n", elementCount, max, min);
        for(uiint i=0; i<epartTable.at(iRegion).size(); i++) {
            if(epartTable.at(iRegion).at(i) == iRegion) {
                uiint id = pMW->getElementID(i);
                pMW->SelectElement_ID(id);
                iint* nodeIDArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                pMW->GetElementVertNodeID(nodeIDArray);
                fprintf(fp, "Tetra %d ", id);
                for(uiint j=0; j<pMW->GetNumOfElementVert(); j++) {
                    fprintf(fp, "%d ", nodeIDArray[j]);
                }
                fprintf(fp, "\n");
            }
        }
        fprintf(fp, "End\n\n");

        //Element group
        fprintf(fp, "ElementGroup\n");
        for(uiint i=0; i<pMW->GetNumOfElementGroup(); i++) {
            string groupName;
            groupName = pMW->GetElementGroupName(i);
            fprintf(fp, "%d %s 0\n", i, groupName.c_str());
        }
        fprintf(fp, "End\n\n");

        //Element group entity
        for(uiint i=0; i<pMW->GetNumOfElementGroup(); i++) {
            fprintf(fp, "ElementGroupEntity\n");
            fprintf(fp, "%d %d", i, 0);
            uiint count = 0;
            for(uiint j=0; j<pMW->GetNumOfElementID(i); j++) {
                uiint element_id = pMW->GetElementID_with_ElementGroup(i, j);
                uiint element_index = pMW->getElementIndex(element_id);
                if(epartTable.at(iRegion).at(element_index) == iRegion) {
                    if(count%10 == 0) {
                        fprintf(fp, "\n");
                    }
                    fprintf(fp, "%d ", element_id);
                    count++;
                }
            }
            fprintf(fp, "\nEnd\n\n");
        }


        //boundary node
        fprintf(fp, "BoundaryNodeMesh\n");
        fprintf(fp, "0 %d\n", pMW->GetNumOfBoundaryNodeMesh());
        for(uiint i=0; i<pMW->GetNumOfBoundaryNodeMesh(); i++) {
            uiint type = pMW->GetBNDType_BNodeMesh(i);
            string name = pMW->GetBNodeMesh_Name(i);
            string typeName;
            if(type == pMW->getNeumannType()) {
                typeName = "Neumann";
            } else {
                typeName = "Dirichlet";
            }
            fprintf(fp, "%d %s %s\n", i, typeName.c_str(), name.c_str());
        }
        fprintf(fp, "End\n\n");
        for(uiint i=0; i<pMW->GetNumOfBoundaryNodeMesh(); i++) {
            fprintf(fp, "BoundaryNode\n");
            uiint type = pMW->GetBNDType_BNodeMesh(i);
            uiint num = pMW->GetNumOfBNode_BNodeMesh(i);
            uiint count = 0;
            string typeName;
            for(uiint j=0; j<num; j++) {
                uiint id = pMW->GetNodeID_BNode_BNodeMesh(i, j);
                uiint gIndex = pMW->getNodeIndex(id);
                uiint region = npartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    count++;
                }
            }
            if(type == pMW->getNeumannType()) {
                typeName = "Neumann";
            } else {
                typeName = "Dirichlet";
            }

            fprintf(fp, "%d %s 0 %d\n", i, typeName.c_str(), count);
            uiint bnodeIndex = 0;
            for(uiint j=0; j<num; j++) {
                uiint id = pMW->GetNodeID_BNode_BNodeMesh(i, j);
                uiint gIndex = pMW->getNodeIndex(id);
                uiint region = npartTable.at(iRegion).at(gIndex);
                double x,y,z;
                pMW->GetNodeCoord(id, x, y, z);
                if(region == iRegion) {
                    for(uiint dof=0; dof<pMW->GetNumOfDOF_BNodeMesh(i, j); dof++) {
                        fprintf(fp, "%d %d %f %f %f ", bnodeIndex, id, x, y, z);
                        fprintf(fp, "%d ", pMW->GetDOF_BNodeMesh(i,j,dof));
                        fprintf(fp, "%f ", pMW->GetBNodeValue_BNodeMesh(i, j, pMW->GetDOF_BNodeMesh(i,j,dof)));
                        fprintf(fp, "\n");
                    }
                    bnodeIndex++;
                }
            }
            fprintf(fp, "End\n\n");
        }

        //boundary volume
        fprintf(fp, "BoundaryVolumeMesh\n");
        fprintf(fp, "0 %d\n", pMW->GetNumOfBoundaryVolumeMesh());
        for(uiint i=0; i<pMW->GetNumOfBoundaryVolumeMesh(); i++) {
            uiint type = pMW->GetBNDType_BVolumeMesh(i);
            string name = pMW->GetBVolumeMesh_Name(i);
            uiint nDOF = pMW->GetNumOfDOF_BVolumeMesh(i);
            string typeName;
            if(type == pMW->getNeumannType()) {
                typeName = "Neumann";
            } else {
                typeName = "Dirichlet";
            }
            fprintf(fp, "%d %s %s %d ", i, typeName.c_str(), name.c_str(), nDOF);
            for(uiint j=0; j<nDOF; j++) {
                fprintf(fp, "%d ", pMW->GetDOF_BVolumeMesh(i, j));
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "End\n\n");


        for(uiint i=0; i<pMW->GetNumOfBoundaryVolumeMesh(); i++) {
            fprintf(fp, "BoundaryVolume\n");
            uiint type = pMW->GetBNDType_BVolumeMesh(i);
            uiint num = pMW->GetNumOfBNode_BVolumeMesh(i);
            uiint count = 0;
            uiint nodecount = 0;
            for(uiint j=0; j<pMW->GetNumOfBVolume(i); j++) {
                uiint id = pMW->GetElemID_BVolume(i, j);
                uiint gIndex = pMW->getElementIndex(id);
                uiint region = epartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    count++;
                }
            }
            for(uiint j=0; j<pMW->GetNumOfBNode_BVolumeMesh(i); j++) {
                uiint id = pMW->GetNodeID_BNode_BVolumeMesh(i, j);
                uiint gIndex = pMW->getNodeIndex(id);
                uiint region = npartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    nodecount++;
                }
            }
            string typeName;
            if(type == pMW->getNeumannType()) {
                typeName = "Neumann";
            } else {
                typeName = "Dirichlet";
            }

            fprintf(fp, "%d %s 0 %d %d\n", i, typeName.c_str(), nodecount, count);
            uiint bvolIndex = 0;
            map<uiint, uiint> bvolNodeTable; //nodeID, bvolIndex;
            for(uiint j=0; j<pMW->GetNumOfBNode_BVolumeMesh(i); j++) {
                uiint id = pMW->GetNodeID_BNode_BVolumeMesh(i, j);
                uiint gIndex = pMW->getNodeIndex(id);
                uiint region = npartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    if(type == pMW->getNeumannType()) {
                        fprintf(fp, "%d %d\n", bvolIndex, id);
                    } else {
                        fprintf(fp, "%d %d ", bvolIndex, id);
                        fprintf(fp, "%d ", pMW->GetNumOfDOF_BVolumeMesh(i));
                        for(uiint dof=0; dof<pMW->GetNumOfDOF_BVolumeMesh(i); dof++) {
                            fprintf(fp, "%d ", pMW->GetDOF_BVolumeMesh(i, dof));
                            fprintf(fp, "%f ", pMW->GetBNodeValue_BVolumeMesh(i, j, pMW->GetDOF_BVolumeMesh(i, j), 0));
                        }
                        fprintf(fp, "\n");
                    }
                    bvolNodeTable.insert( map<uiint, uiint>::value_type(id, bvolIndex) );
                    bvolIndex++;
                }
            }
            bvolIndex = 0;
            for(uiint j=0; j<pMW->GetNumOfBVolume(i); j++) {
                uiint id = pMW->GetElemID_BVolume(i, j);
                uiint gIndex = pMW->getElementIndex(id);
                uiint region = epartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    if(type == pMW->getNeumannType()) {
                        for(uiint dof=0; dof<pMW->GetNumOfDOF_BVolumeMesh(i); dof++) {
                            fprintf(fp, "Tetra %d %d 0 ", bvolIndex, id);
                            fprintf(fp, "%d ", pMW->GetDOF_BVolumeMesh(i, dof));
                            pMW->SelectElement_ID(id);
                            iint* nodeIDArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                            pMW->GetElementVertNodeID(nodeIDArray);
                            for(uiint k=0; k<pMW->GetNumOfElementVert(); k++) {
                                fprintf(fp, "%d ", bvolNodeTable[nodeIDArray[k]]);
                            }
                            free(nodeIDArray);
                            fprintf(fp, "%f ", pMW->GetBVolumeValue(i, j, pMW->GetDOF_BVolumeMesh(i, j)));
                            fprintf(fp, "\n");
                        }
                        bvolIndex++;
                    } else {
                        fprintf(fp, "Tetra %d %d 0 ", bvolIndex, id);
                        pMW->SelectElement_ID(id);
                        iint* nodeIDArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                        pMW->GetElementVertNodeID(nodeIDArray);
                        for(uiint k=0; k<pMW->GetNumOfElementVert(); k++) {
                            fprintf(fp, "%d ", bvolNodeTable[nodeIDArray[k]]);
                        }
                        free(nodeIDArray);
                        fprintf(fp, "\n");
                        bvolIndex++;
                    }
                }
            }
            fprintf(fp, "End\n\n");
        }


        //boundary face
        fprintf(fp, "BoundaryFaceMesh\n");
        fprintf(fp, "0 %d\n", pMW->GetNumOfBoundaryFaceMesh());
        for(uiint i=0; i<pMW->GetNumOfBoundaryFaceMesh(); i++) {
            uiint type = pMW->GetBNDType_BFaceMesh(i);
            string name = pMW->GetBFaceMesh_Name(i);
            uiint nDOF = pMW->GetNumOfDOF_BFaceMesh(i);
            string typeName;
            if(type == pMW->getNeumannType()) {
                typeName = "Neumann";
            } else {
                typeName = "Dirichlet";
            }
            fprintf(fp, "%d %s %s %d ", i, typeName.c_str(), name.c_str(), nDOF);
            for(uiint j=0; j<nDOF; j++) {
                fprintf(fp, "%d ", pMW->GetDOF_BFaceMesh(i, j));
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "End\n\n");


        for(uiint i=0; i<pMW->GetNumOfBoundaryFaceMesh(); i++) {
            fprintf(fp, "BoundaryFace\n");
            uiint type = pMW->GetBNDType_BFaceMesh(i);
            uiint num = pMW->GetNumOfBNode_BFaceMesh(i);
            uiint count = 0;
            uiint nodecount = 0;
            for(uiint j=0; j<pMW->GetNumOfBFace(i); j++) {
                uiint id = pMW->GetElemID_BFace(i, j);
                uiint gIndex = pMW->getElementIndex(id);
                uiint region = epartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    count++;
                }
            }
            for(uiint j=0; j<pMW->GetNumOfBNode_BFaceMesh(i); j++) {
                uiint id = pMW->GetNodeID_BNode_BFaceMesh(i, j);
                uiint gIndex = pMW->getNodeIndex(id);
                uiint region = npartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    nodecount++;
                }
            }
            string typeName;
            if(type == pMW->getNeumannType()) {
                typeName = "Neumann";
            } else {
                typeName = "Dirichlet";
            }

            fprintf(fp, "%d %s 0 %d %d\n", i, typeName.c_str(), nodecount, count);
            uiint bfaceIndex = 0;
            map<uiint, uiint> bfaceNodeTable; //nodeID, bfaceIndex;
            for(uiint j=0; j<pMW->GetNumOfBNode_BFaceMesh(i); j++) {
                uiint id = pMW->GetNodeID_BNode_BFaceMesh(i, j);
                uiint gIndex = pMW->getNodeIndex(id);
                uiint region = npartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    fprintf(fp, "%d %d ", bfaceIndex, id);
                    if(type == pMW->getNeumannType()) {

                    } else {
                        fprintf(fp, "%d ", pMW->GetNumOfDOF_BFaceMesh(i));
                        for(uiint dof=0; dof<pMW->GetNumOfDOF_BFaceMesh(i); dof++) {
                            fprintf(fp, "%d ", pMW->GetDOF_BFaceMesh(i, dof));
                            fprintf(fp, "%f ", pMW->GetBNodeValue_BFaceMesh(i, j, pMW->GetDOF_BFaceMesh(i, dof), 0));
                        }
                    }
                    fprintf(fp, "\n");
                    bfaceNodeTable.insert( map<uiint, uiint>::value_type(id, bfaceIndex) );
                    bfaceIndex++;
                }
            }
            bfaceIndex = 0;
            for(uiint j=0; j<pMW->GetNumOfBFace(i); j++) {
                uiint id = pMW->GetElemID_BFace(i, j);
                uiint gIndex = pMW->getElementIndex(id);
                uiint region = epartTable.at(iRegion).at(gIndex);
                uiint faceId = pMW->GetFaceID_BFace(i, j);
                if(region == iRegion) {
                    if(type == pMW->getNeumannType()) {
                        for(uiint dof=0; dof<pMW->GetNumOfDOF_BFaceMesh(i); dof++) {
                            fprintf(fp, "Triangle %d %d %d ", bfaceIndex, id, faceId);
                            fprintf(fp, "%d ", pMW->GetDOF_BFaceMesh(i, dof));
                            pMW->SelectElement_ID(id);
                            iint* nodeIDArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                            pMW->GetElementVertNodeID(nodeIDArray);
                            switch(faceId) {
                            case 0:
                                fprintf(fp, "%d %d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[1]], bfaceNodeTable[nodeIDArray[2]]);
                                break;
                            case 1:
                                fprintf(fp, "%d %d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[3]], bfaceNodeTable[nodeIDArray[1]]);
                                break;
                            case 2:
                                fprintf(fp, "%d %d %d ", bfaceNodeTable[nodeIDArray[1]], bfaceNodeTable[nodeIDArray[3]], bfaceNodeTable[nodeIDArray[2]]);
                                break;
                            case 3:
                                fprintf(fp, "%d %d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[2]], bfaceNodeTable[nodeIDArray[3]]);
                                break;
                            }
                            free(nodeIDArray);
                            fprintf(fp, "%f ", pMW->GetBFaceValue(i, j, pMW->GetDOF_BFaceMesh(i, dof)));
                            fprintf(fp, "\n");
                        }
                        bfaceIndex++;
                    } else {
                        fprintf(fp, "Triangle %d %d %d ", bfaceIndex, id, faceId);
                        pMW->SelectElement_ID(id);
                        iint* nodeIDArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                        pMW->GetElementVertNodeID(nodeIDArray);
                        switch(faceId) {
                        case 0:
                            fprintf(fp, "%d %d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[1]], bfaceNodeTable[nodeIDArray[2]]);
                            break;
                        case 1:
                            fprintf(fp, "%d %d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[3]], bfaceNodeTable[nodeIDArray[1]]);
                            break;
                        case 2:
                            fprintf(fp, "%d %d %d ", bfaceNodeTable[nodeIDArray[1]], bfaceNodeTable[nodeIDArray[3]], bfaceNodeTable[nodeIDArray[2]]);
                            break;
                        case 3:
                            fprintf(fp, "%d %d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[2]], bfaceNodeTable[nodeIDArray[3]]);
                            break;
                        }
                        free(nodeIDArray);
                        fprintf(fp, "\n");
                        bfaceIndex++;
                    }


                }
            }
            fprintf(fp, "End\n\n");
        }

        //boundary edge
        fprintf(fp, "BoundaryEdgeMesh\n");
        fprintf(fp, "0 %d\n", pMW->GetNumOfBoundaryEdgeMesh());
        for(uiint i=0; i<pMW->GetNumOfBoundaryEdgeMesh(); i++) {
            uiint type = pMW->GetBNDType_BEdgeMesh(i);
            string name = pMW->GetBEdgeMesh_Name(i);
            uiint nDOF = pMW->GetNumOfDOF_BEdgeMesh(i);
            string typeName;
            if(type == pMW->getNeumannType()) {
                typeName = "Neumann";
            } else {
                typeName = "Dirichlet";
            }
            fprintf(fp, "%d %s %s %d ", i, typeName.c_str(), name.c_str(), nDOF);
            for(uiint j=0; j<nDOF; j++) {
                fprintf(fp, "%d ", pMW->GetDOF_BEdgeMesh(i, j));
            }
            fprintf(fp, "0\n");
        }
        fprintf(fp, "End\n\n");


        for(uiint i=0; i<pMW->GetNumOfBoundaryEdgeMesh(); i++) {
            fprintf(fp, "BoundaryEdge\n");
            uiint type = pMW->GetBNDType_BEdgeMesh(i);
            uiint num = pMW->GetNumOfBNode_BEdgeMesh(i);
            uiint count = 0;
            uiint nodecount = 0;
            for(uiint j=0; j<pMW->GetNumOfBEdge(i); j++) {
                uiint id = pMW->GetElemID_BEdge(i, j);
                uiint gIndex = pMW->getElementIndex(id);
                uiint region = epartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    count++;
                }
            }
            for(uiint j=0; j<pMW->GetNumOfBNode_BEdgeMesh(i); j++) {
                uiint id = pMW->GetNodeID_BNode_BEdgeMesh(i, j);
                uiint gIndex = pMW->getNodeIndex(id);
                uiint region = npartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    nodecount++;
                }
            }
            string typeName;
            if(type == pMW->getNeumannType()) {
                typeName = "Neumann";
            } else {
                typeName = "Dirichlet";
            }

            fprintf(fp, "%d %s 0 %d %d\n", i, typeName.c_str(), nodecount, count);
            uiint bfaceIndex = 0;
            map<uiint, uiint> bfaceNodeTable; //nodeID, bfaceIndex;
            for(uiint j=0; j<pMW->GetNumOfBNode_BEdgeMesh(i); j++) {
                uiint id = pMW->GetNodeID_BNode_BEdgeMesh(i, j);
                uiint gIndex = pMW->getNodeIndex(id);
                uiint region = npartTable.at(iRegion).at(gIndex);
                if(region == iRegion) {
                    fprintf(fp, "%d %d ", bfaceIndex, id);
                    if(type == pMW->getNeumannType()) {
                    } else {
                        fprintf(fp, "%d ", pMW->GetNumOfDOF_BEdgeMesh(i));
                        for(uiint dof=0; dof<pMW->GetNumOfDOF_BEdgeMesh(i); dof++) {
                            fprintf(fp, "%d ", pMW->GetDOF_BEdgeMesh(i, dof));
                            fprintf(fp, "%f ", pMW->GetBNodeValue_BEdgeMesh(i, j, pMW->GetDOF_BEdgeMesh(i, j), 0));
                        }
                    }
                    fprintf(fp, "\n");
                    bfaceNodeTable.insert( map<uiint, uiint>::value_type(id, bfaceIndex) );
                    bfaceIndex++;
                }
            }
            bfaceIndex = 0;
            for(uiint j=0; j<pMW->GetNumOfBEdge(i); j++) {
                uiint id = pMW->GetElemID_BEdge(i, j);
                uiint gIndex = pMW->getElementIndex(id);
                uiint region = epartTable.at(iRegion).at(gIndex);
                uiint faceId = pMW->GetEdgeID_BEdge(i, j);
                if(region == iRegion) {
                    if(type == pMW->getNeumannType()) {
                        for(uiint dof=0; dof<pMW->GetNumOfDOF_BEdgeMesh(i); dof++) {
                            fprintf(fp, "Beam %d %d %d ", bfaceIndex, id, faceId);
                            fprintf(fp, "%d ", pMW->GetDOF_BEdgeMesh(i, dof));
                            pMW->SelectElement_ID(id);
                            iint* nodeIDArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                            pMW->GetElementVertNodeID(nodeIDArray);
                            switch(faceId) {
                            case 0:
                                fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[1]]);
                                break;
                            case 1:
                                fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[1]], bfaceNodeTable[nodeIDArray[2]]);
                                break;
                            case 2:
                                fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[2]], bfaceNodeTable[nodeIDArray[0]]);
                                break;
                            case 3:
                                fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[3]]);
                                break;
                            case 4:
                                fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[1]], bfaceNodeTable[nodeIDArray[3]]);
                                break;
                            case 5:
                                fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[2]], bfaceNodeTable[nodeIDArray[3]]);
                                break;
                            }
                            free(nodeIDArray);
                            fprintf(fp, "%f ", pMW->GetBEdgeValue(i, j, pMW->GetDOF_BEdgeMesh(i, j)));
                            fprintf(fp, "\n");
                        }
                        bfaceIndex++;
                    } else {
                        fprintf(fp, "Beam %d %d %d ", bfaceIndex, id, faceId);
                        pMW->SelectElement_ID(id);
                        iint* nodeIDArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                        pMW->GetElementVertNodeID(nodeIDArray);
                        switch(faceId) {
                        case 0:
                            fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[1]]);
                            break;
                        case 1:
                            fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[1]], bfaceNodeTable[nodeIDArray[2]]);
                            break;
                        case 2:
                            fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[2]], bfaceNodeTable[nodeIDArray[0]]);
                            break;
                        case 3:
                            fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[0]], bfaceNodeTable[nodeIDArray[3]]);
                            break;
                        case 4:
                            fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[1]], bfaceNodeTable[nodeIDArray[3]]);
                            break;
                        case 5:
                            fprintf(fp, "%d %d ", bfaceNodeTable[nodeIDArray[2]], bfaceNodeTable[nodeIDArray[3]]);
                            break;
                        }
                        free(nodeIDArray);
                        fprintf(fp, "\n");
                        bfaceIndex++;
                    }

                }
            }
            fprintf(fp, "End\n\n");
        }


        fprintf(fp, "CommMesh2\n");
        int num = 0;
        for(int i=0; i<nparts; i++) {
            if(commFaceNumTable.at(iRegion).at(i) > 0) {
                num++;
            }
        }
        fprintf(fp, "0 %d\n", num);
        for(int i=0; i<nparts; i++) {
            if(commFaceNumTable.at(iRegion).at(i) > 0) {
                int numNode = 0;
                for(int j=0; j<commNodeTable.at(iRegion).at(i).size(); j++) {
                    if(commNodeTable.at(iRegion).at(i).at(j) == i) {
                        numNode++;
                    }
                }
                //fprintf(fp, "%d %d %d %d %d\n", getCommMeshID(iRegion, i, nparts), commFaceNum.at(i), numNode, iRegion, i);
                fprintf(fp, "%d %d %d %d %d\n", this->getCommMeshIX(iRegion, i), commFaceNumTable.at(iRegion).at(i), numNode, iRegion, i);
            }
        }
        fprintf(fp, "End\n\n");

        for(int i=0; i<nparts; i++) {
            if(commFaceNumTable.at(iRegion).at(i) > 0) {
                fprintf(fp, "CommNodeCM2\n");
                int numNode = 0;
                for(int j=0; j<commNodeTable.at(iRegion).at(i).size(); j++) {
                    if(commNodeTable.at(iRegion).at(i).at(j) == i) {
                        numNode++;
                    }
                }
                fprintf(fp, "%d 0 %d\n", this->getCommMeshIX(iRegion, i), numNode);
                for(int j=0; j<commNodeTable.at(iRegion).at(i).size(); j++) {
                    if(commNodeTable.at(iRegion).at(i).at(j) == i) {
                        double x,y,z;
                        pMW->GetNodeCoord(pMW->getNodeID(j), x, y, z);
                        fprintf(fp, "%d %d %f %f %f\n", commNodeIndexTable.at(iRegion).at(i).at(j), pMW->getNodeID(j), x, y, z);
                    }
                }
                fprintf(fp, "End\n\n");
            }
        }

        for(int i=0; i<nparts; i++) {
            if(commFaceNumTable.at(iRegion).at(i) > 0) {
                fprintf(fp, "CommFace\n");
                fprintf(fp, "%d 0 %d\n", this->getCommMeshIX(iRegion, i), commFaceNumTable.at(iRegion).at(i));
                int commFaceI = 0;
                for(int j=0; j<commFaceRegionTable.at(iRegion).size(); j++) {
                    if(commFaceRegionTable.at(iRegion).at(j) == i) {
                        uiint id,face,type;
                        iint* nodeIXArray;

                        if(iRegion < i) {
                            id = pMW->getElementID(commFaceElementTable.at(iRegion).at(j));
                            face = commFaceIndexTable.at(iRegion).at(j);
                            type = commFaceTypeTable.at(iRegion).at(j);
                            pMW->SelectElement_ID(id);
                            nodeIXArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                            pMW->GetElementVertNodeIndex(nodeIXArray);
                        } else {
                            id = pMW->getElementID(commFaceElementTable.at(iRegion).at(j));
                            face = commFaceIndexVertTable.at(iRegion).at(j);
                            type = commFaceTypeTable.at(iRegion).at(j);
                            pMW->SelectElement_ID(pMW->getElementID(commFaceElementVertTable.at(iRegion).at(j)));
                            nodeIXArray = (iint*)malloc(sizeof(iint) * pMW->GetNumOfElementVert());
                            pMW->GetElementVertNodeIndex(nodeIXArray);
                        }

                        if(type == 0) {
                            fprintf(fp, "Point %d %d %d ",  commFaceI, id, commFaceIndexTable.at(iRegion).at(j));
                            fprintf(fp, "%d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[face]));
                        } else if(type == 1) {
                            fprintf(fp, "Beam %d %d %d ",  commFaceI, id, commFaceIndexTable.at(iRegion).at(j));
                            switch(face) {
                            case 0:
                                fprintf(fp, "%d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[0]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[1]));
                                break;
                            case 1:
                                fprintf(fp, "%d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[1]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[2]));
                                break;
                            case 2:
                                fprintf(fp, "%d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[2]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[0]));
                                break;
                            case 3:
                                fprintf(fp, "%d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[0]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[3]));
                                break;
                            case 4:
                                fprintf(fp, "%d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[1]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[3]));
                                break;
                            case 5:
                                fprintf(fp, "%d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[2]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[3]));
                                break;
                            }
                        } else if(type == 3) {
                            fprintf(fp, "Triangle %d %d %d ",  commFaceI, id, commFaceIndexTable.at(iRegion).at(j));
                            switch(face) {
                            case 0:
                                fprintf(fp, "%d %d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[0]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[1]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[2]));
                                break;
                            case 1:
                                fprintf(fp, "%d %d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[0]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[3]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[1]));
                                break;
                            case 2:
                                fprintf(fp, "%d %d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[1]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[3]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[2]));
                                break;
                            case 3:
                                fprintf(fp, "%d %d %d\n", commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[0]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[2]), commNodeIndexTable.at(iRegion).at(i).at(nodeIXArray[3]));
                                break;
                            }
                        }
                        commFaceI++;
                        free(nodeIXArray);
                    }
                }
                fprintf(fp, "End\n\n");
            }
        }



        fclose(fp);
    }

}

void CPartitioner::writeUCD(string filename)
{
}

void CPartitioner::printReport()
{
    printf("\n");
    printf("--------------------------------------------------------------------------------\n");
    printf("Input mesh\n");
    printf("--------------------------------------------------------------------------------\n");
    printf("number of nodes    : %d\n", nn);
    printf("number of elenents : %d\n\n\n", ne);

    for(uiint iRegion=0; iRegion < nparts; iRegion++) {
        printf("--------------------------------------------------------------------------------\n");
        printf("Rank #%d\n", iRegion);
        printf("--------------------------------------------------------------------------------\n");
        printf("number of nodes    : %d\n", numNodesTable.at(iRegion));
        printf("number of elenents : %d\n", numElementsTable.at(iRegion));
        int numNode = 0;
        for(int i=0; i<nparts; i++) {
            if(commFaceNumTable.at(iRegion).at(i) > 0) {
                for(int j=0; j<commNodeTable.at(iRegion).at(i).size(); j++) {
                    if(commNodeTable.at(iRegion).at(i).at(j) == i) {
                        numNode++;
                    }
                }
            }
        }
        printf("number of commNodes: %d\n\n\n", numNode);
    }

}

uiint CPartitioner::getCommMeshIX(uiint myrank, uiint rank)
{
    uiint x,y;
    if(myrank < rank) {
        x = rank;
        y = myrank;
    } else if(rank < myrank) {
        x = myrank;
        y = rank;
    } else {
        return -1;;
    }

    uiint id = -1;
    uiint i,j;
    for(i=0; i<nparts; i++) {
        for(j=0; j<nparts; j++) {
            if(j > i) {
                if(globalComm.at(getCommMeshID(i, j, nparts)) > 0) {
                    id++;
                    if(i == y && j == x) {
                        return id;
                    }
                } else {

                }
            }
        }
    }
    return -1;
}
