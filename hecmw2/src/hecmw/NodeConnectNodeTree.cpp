//
//  NodeConnectNodeTree.cpp
//
//
//
//                  2009.08.11
//                  2009.08.11
//                  k.Takeda
#include "NodeConnectNodeTree.h"
using namespace pmw;

CNodeConnectNodeTree::CNodeConnectNodeTree()
{
    // Hexa
    // --
    uint hexcon[8][3]={
    {1,3,4},
    {0,5,2},
    {1,3,6},
    {0,2,7},
    {0,5,7},
    {4,1,6},
    {5,2,7},
    {4,3,6}
    };
    
    uint i,ii;
    mvHexaConnectNode.resize(8);
    for(i=0; i< 8; i++){
        mvHexaConnectNode[i].resize(3);
    };
    for(i=0; i< 8; i++){
        for(ii=0; ii< 3; ii++){
            mvHexaConnectNode[i][ii]= hexcon[i][ii];
        };
    };
    ////debug
    //cout << "CNodeConnectNodeTree() Hexa" << endl;


    // Tetra
    // --
    uint tetcon[4][3]={
    {2,3,1},
    {0,3,2},
    {1,3,0},
    {2,0,1}
    };
    
    mvTetraConnectNode.resize(4);
    for(i=0; i< 4; i++){
        mvTetraConnectNode[i].resize(3);
    };
    for(i=0; i< 4; i++){
        for(ii=0; ii< 3; ii++){
            mvTetraConnectNode[i][ii]= tetcon[i][ii];
        };
    };
    ////debug
    //cout << "CNodeConnectNodeTree() Tetra" << endl;


    // Prism
    // --
    uint pricon[6][3]={
    {2,3,1},
    {0,4,2},
    {1,5,0},
    {5,0,4},
    {3,1,5},
    {4,3,2}
    };
    
    mvPrismConnectNode.resize(6);
    for(i=0; i< 6; i++){
        mvPrismConnectNode[i].resize(3);
    };
    for(i=0; i< 6; i++){
        for(ii=0; ii< 3; ii++){
            mvPrismConnectNode[i][ii]= pricon[i][ii];
        };
    };
    ////debug
    //cout << "CNodeConnectNodeTree() Prism" << endl;

    
    // Pyramid
    // --
    uint pycon[5][4]={
    {1,4,3,1},// 本来は3つまで. 最後の値は先頭と同値
    {0,4,2,0},//  ↑ 同上
    {1,4,3,1},//  ↑ 同上
    {0,4,2,0},//  ↑ 同上
    {0,1,2,3} // ４節点と接続
    };
    
    mvPyramidConnectNode.resize(5);
    for(i=0; i< 5; i++){
        if(i!=4){
            mvPyramidConnectNode[i].resize(3);
        }else{
            mvPyramidConnectNode[i].resize(4);
        }
    };
    for(i=0; i< 5; i++){
        if(i!=4) for(ii=0; ii< 3; ii++) mvPyramidConnectNode[i][ii]= pycon[i][ii];
        if(i==4) for(ii=0; ii< 4; ii++) mvPyramidConnectNode[i][ii]= pycon[i][ii];
    };
    ////debug
    //cout << "CNodeConnectNodeTree() Pyramid" << endl;

    
    // Quad
    // --
    uint qucon[4][2]={
    {1,3},
    {0,2},
    {1,3},
    {2,0}
    };

    mvQuadConnectNode.resize(4);
    for(i=0; i< 4; i++){
        mvQuadConnectNode[i].resize(2);
    };
    for(i=0; i< 4; i++){
        for(ii=0; ii< 2; ii++){
            mvQuadConnectNode[i][ii]= qucon[i][ii];
        };
    };
    ////debug
    //cout << "CNodeConnectNodeTree() Quad" << endl;


    // Triangle
    // --
    uint tricon[3][2]={
        {1,2},
        {0,2},
        {0,1}
    };

    mvTriangleConnectNode.resize(3);
    for(i=0; i< 3; i++){
        mvTriangleConnectNode[i].resize(2);
    };
    for(i=0; i< 3; i++){
        for(ii=0; ii< 2; ii++){
            mvTriangleConnectNode[i][ii]= tricon[i][ii];
        };
    };
    ////debug
    //cout << "CNodeConnectNodeTree() Triangle" << endl;

    // Beam
    // --
    uint beamcon[2][1]={
        {1},
        {0}
    };

    mvBeamConnectNode.resize(2);
    for(i=0; i< 2; i++){
        mvBeamConnectNode[i].resize(1);
    };
    for(i=0; i< 2; i++){
        for(ii=0; ii< 1; ii++){
            mvBeamConnectNode[i][ii]= beamcon[i][ii];
        };
    };
    ////debug
    //cout << "CNodeConnectNodeTree() Beam" << endl;

}

CNodeConnectNodeTree::~CNodeConnectNodeTree()
{
    ;
}




