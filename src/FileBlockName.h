//
//  FileBlockName  struct
// 
//
//				2009.09.18
//				2008.12.10
//				k.Takeda
#ifndef FILE_BLOCK_NAME_7B85EDDD_CCF3_4c0b_B435_154AA1897A9C
#define FILE_BLOCK_NAME_7B85EDDD_CCF3_4c0b_B435_154AA1897A9C

/*
	<< HEC_MW2 Tag >>
	<<  FileTagName.h  >>

	<< cnt file  >>
	!CONTROL
	!MESH
	!MESH GROUP
	!RESTART
	!RESULT

	<< mesh file >>
	!HEADER
	!ZERO
	!NODE
	!ELEMENT
	!NGROUP
	!EGROUP
	!SGROUP
        !LGROUP //新グループ
	!EQUATION
	!AMPLITUDE
	!SECTION
	!MATERIAL
	!INITIAL CONDITION
	!INCLUDE
	!CONNECTIVITY
	!END
*/

namespace FileIO{
struct FileBlockName{
    // Cnt
    // ----
    static const char* StartMeshFileName(){ return "MeshFileName";}
    static const char* EndMeshFileName(){ return "End";}


    // Mesh
    // ----
    // Node
    static const char* StartNode(){ return "Node";}
    static const char* EndNode(){ return "End";}
    // Element
    static const char* StartElement(){ return "Element";}
    static const char* EndElement(){ return "End";}
    
    
    // Boundary
    // ----
    // BounaryNodeMesh(節点型境界の境界ID全体数)
    static const char* StartBoundaryNodeMesh(){ return "BoundaryNodeMesh";}
    static const char* EndBoundaryNodeMesh(){ return "End";}
    // BounaryFaceMesh(面型境界の境界ID全体数)
    static const char* StartBoundaryFaceMesh(){ return "BoundaryFaceMesh";}
    static const char* EndBoundaryFaceMesh(){ return "End";}
    // BounaryVolumeMesh(体積型境界の境界ID全体数)
    static const char* StartBoundaryVolumeMesh(){ return "BoundaryVolumeMesh";}
    static const char* EndBoundaryVolumeMesh(){ return "End";}
    // BounaryEdgeMesh(辺型境界の境界ID全体数)
    static const char* StartBoundaryEdgeMesh(){ return "BoundaryEdgeMesh";}
    static const char* EndBoundaryEdgeMesh(){ return "End";}

    // BoundaryNode
    static const char* StartBoundaryNode(){ return "BoundaryNode";}
    static const char* EndBoundaryNode(){ return "End";}
    // BoundaryFace
    static const char* StartBoundaryFace(){ return "BoundaryFace";}
    static const char* EndBoundaryFace(){ return "End";}
    // BounaryEdge
    static const char* StartBoundaryEdge(){ return "BoundaryEdge";}
    static const char* EndBoundaryEdge(){ return "End";}
    // BoundaryVolume
    static const char* StartBoundaryVolume(){ return "BoundaryVolume";}
    static const char* EndBoundaryVolume(){ return "End";}
    


    // Material
    // ----
    static const char* StartMaterial(){ return "Material";}
    static const char* EndMaterial(){ return "End";}


    // Group
    // ----
    static const char* StartElementGroup(){ return "ElementGroup";}
    static const char* EndElementGroup(){ return "End";}

    static const char* StartElementGroupEntity(){ return "ElementGroupEntity";}
    static const char* EndElementGroupEntity(){ return "End";}



    // Communication Type 1
    // ----
    // CommMesh(節点分割型==要素共有型)
    // ----
    static const char* StartCommMesh(){ return "CommMesh";}
    static const char* EndCommMesh(){ return "End";}
    // CommNode
    static const char* StartCommNode(){ return "CommNode";}
    static const char* EndCommNode(){ return "End";}
    // CommElement
    static const char* StartCommElement(){ return "CommElement";}
    static const char* EndCommElement(){ return "End";}

    // Communication Type 2
    // ----
    // CommMesh2 (要素分割型==節点共有型)
    // ----
    static const char* StartCommMesh2(){ return "CommMesh2";}
    static const char* EndCommMesh2(){ return  "End";}
    // CommNode for CommMesh2
    static const char* StartCommNodeCM2(){ return "CommNodeCM2";}
    static const char* EndCommNodeCM2(){ return "End";}
    // CommFace
    static const char* StartCommFace(){ return "CommFace";}
    static const char* EndCommFace(){ return "End";}




    // GMGModel(MultiGrid Model)
    // ----
    static const char* StartRefine(){ return "Refine";}
    static const char* EndRefine(){ return "End";}

    
    // AssyModel(Assembly Model)
    // ----
    static const char* StartAssyModel(){ return "AssyModel";}
    static const char* EndAssyModel(){ return "End";}
    // Mesh(Part)
    static const char* StartMesh(){ return "Mesh";}
    static const char* EndMesh(){ return "End";}
    // ContactMesh(Interface)
    static const char* StartContactMesh(){ return "ContactMesh";}
    static const char* EndContactMesh(){ return "End";}
    
};
}
#endif
