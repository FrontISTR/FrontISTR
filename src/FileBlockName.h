/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileBlockName.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef FILE_BLOCK_NAME_7B85EDDD_CCF3_4c0b_B435_154AA1897A9C
#define FILE_BLOCK_NAME_7B85EDDD_CCF3_4c0b_B435_154AA1897A9C
namespace FileIO{
struct FileBlockName{
    static const char* StartMeshFileName(){ return "MeshFileName";}
    static const char* EndMeshFileName(){ return "End";}
    static const char* StartNode(){ return "Node";}
    static const char* EndNode(){ return "End";}
    static const char* StartElement(){ return "Element";}
    static const char* EndElement(){ return "End";}
    static const char* StartBoundaryNode(){ return "BoundaryNode";}
    static const char* EndBoundaryNode(){ return "End";}
    static const char* StartBoundaryFace(){ return "BoundaryFace";}
    static const char* EndBoundaryFace(){ return "End";}
    static const char* StartBoundaryVolume(){ return "BoundaryVolume";}
    static const char* EndBoundaryVolume(){ return "End";}
    static const char* StartMaterial(){ return "Material";}
    static const char* EndMaterial(){ return "End";}
    static const char* StartCommMesh(){ return "CommMesh";}
    static const char* EndCommMesh(){ return "End";}
    static const char* StartCommNode(){ return "CommNode";}
    static const char* EndCommNode(){ return "End";}
    static const char* StartCommElement(){ return "CommElement";}
    static const char* EndCommElement(){ return "End";}
    static const char* StartCommMesh2(){ return "CommMesh2";}
    static const char* EndCommMesh2(){ return  "End";}
    static const char* StartCommNodeCM2(){ return "CommNodeCM2";}
    static const char* EndCommNodeCM2(){ return "End";}
    static const char* StartCommFace(){ return "CommFace";}
    static const char* EndCommFace(){ return "End";}
    static const char* StartRefine(){ return "Refine";}
    static const char* EndRefine(){ return "End";}
    static const char* StartAssyModel(){ return "AssyModel";}
    static const char* EndAssyModel(){ return "End";}
    static const char* StartMesh(){ return "Mesh";}
    static const char* EndMesh(){ return "End";}
    static const char* StartContactMesh(){ return "ContactMesh";}
    static const char* EndContactMesh(){ return "End";}
};
}
#endif
