/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileBlockNameMW2.h
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
#ifndef FILE_TAG_NAME_CC0B8D86_67BD_493a_B552_FBF00C1F9B02
#define FILE_TAG_NAME_CC0B8D86_67BD_493a_B552_FBF00C1F9B02
namespace FileIO
{
struct FileBlockNameMW2 {
    static char HashMark() {
        return '#';
    }
    static char Exclamation() {
        return '!';
    }
    static const char* Header() {
        return "!HEADER";
    }
    static const char* Zero() {
        return "!ZERO";
    }
    static const char* Node() {
        return "!NODE";
    }
    static const char* Element() {
        return "!ELEMENT";
    }
    static const char* NGroup() {
        return "!NGROUP";
    }
    static const char* EGroup() {
        return "!EGROUP";
    }
    static const char* SGroup() {
        return "!SGROUP";
    }
    static const char* Equation() {
        return "!EQUATION";
    }
    static const char* Amplitude() {
        return "!AMPLITUDE";
    }
    static const char* Section() {
        return "!SECTION";
    }
    static const char* Material() {
        return "!MATERIAL";
    }
    static const char* Initial() {
        return "!INITIAL CONDITION";
    }
    static const char* Include() {
        return "!INCLUDE";
    }
    static const char* Connectivity() {
        return "!CONNECTIVITY";
    }
    static const char* End() {
        return "!END";
    }
    static const char* Control() {
        return "!CONTROL";
    }
    static const char* Mesh() {
        return "!MESH";
    }
    static const char* MeshGroup() {
        return "!MESH GROUP";
    }
    static const char* Restart() {
        return "!RESTART";
    }
    static const char* Result() {
        return "!RESULT";
    }
    static const char* Name_Mesh() {
        return "fstrMSH";
    }
    static const char* Name_Control() {
        return "fstrCNT";
    }
    static const char* Name_Result() {
        return "fstrRES";
    }
    static const char* Name_VisMesh() {
        return "mesh";
    }
    static const char* Name_VisIn() {
        return "result";
    }
    static const char* Name_VisOut() {
        return "vis_out";
    }
    static const char* Name_Restart() {
        return "restart";
    }
    static const char* Name_Part_IN() {
        return "part_in";
    }
    static const char* Name_Part_OUT() {
        return "part_out";
    }
    static const char* Type_Dist() {
        return "HECMW-DIST";
    }
    static const char* Type_Entire() {
        return "HECMW-ENTIRE";
    }
    static const char* IO_Out() {
        return "OUT";
    }
    static const char* IO_In() {
        return  "IN";
    }
    static const char* IO_InOut() {
        return "INOUT";
    }
    static const char* Refine() {
        return "!REFINE";
    }
};
struct FileTypeMW2 {
    enum {
        Mesh,
        Control,
        Restart,
        Result,
        PartIN,
        PartOUT,
        VisMesh,
        VisIN,
        VisOUT,
        CADFitFile,
        Limit
    };
};
}
#endif
