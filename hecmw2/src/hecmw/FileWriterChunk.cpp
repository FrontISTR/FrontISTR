/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterChunk.cpp
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
#include "FileWriterChunk.h"
#include "FileBlockName.h"
using namespace FileIO;
CFileWriterChunk::CFileWriterChunk()
{
    mvWriter.reserve(9);
    mvWriter.push_back(new CFileWriterNode);
    mvWriter.push_back(new CFileWriterElement);
    mvWriter.push_back(new CFileWriterBoundaryNode);
    mvWriter.push_back(new CFileWriterBoundaryFace);
    mvWriter.push_back(new CFileWriterBoundaryVolume);
    mvWriter.push_back(new CFileWriterBoundaryEdge);
    mvWriter.push_back(new CFileWriterContactMesh);
    mvWriter.push_back(new CFileWriterCommMesh2);
    mvWriter.push_back(new CFileWriterElementGroup);
    mb_fstr= false;
}
CFileWriterChunk::~CFileWriterChunk()
{
    for_each(mvWriter.begin(), mvWriter.end(), pmw::DeleteObject());
}
void CFileWriterChunk::setSolutionType(const uiint& nSolutionType)
{
    CFileWriter *pWriter;
    uiint numOfWriter=mvWriter.size();
    uiint i;
    for(i=0; i < numOfWriter; i++) {
        pWriter = mvWriter[i];
        pWriter->setSolutionType(nSolutionType);
    }
}

void CFileWriterChunk::WriteDebug(string& filename, const uiint& nNumOfLevel)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::MWDebug,"Debug output filename => ",filename);

    ofstream ofs(filename.c_str(), ios::out);
    uiint i, ilevel;

    for(ilevel=0; ilevel< nNumOfLevel; ilevel++) {
        pLogger->Info(Utility::LoggerMode::Info, "Write, ilevel ", ilevel);

        for(i=0; i < mvWriter.size(); i++) {
            mvWriter[i]->WriteDebug(ofs, ilevel);
        };
    };
    ofs.close();
}

void CFileWriterChunk::markingFstrStyle()
{
    mb_fstr=true;
}
void CFileWriterChunk::WriteRes(const uiint& nStep, string& filename, bool bBinary)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::MWDebug,"Res Filename => ",filename);
    stringstream ss;
    ss << nStep;
    string sFileName;
    if(mb_fstr) {
        sFileName= filename + "." + ss.str();
    } else {
        sFileName= filename + "." + ss.str() + ".res";
    }
    ofstream ofs;
    CFileWriterRes oRes;
    if(bBinary) {
        ofs.open(sFileName.c_str(), ios::out|ios::binary);
        ofs.write(FileBlockName::StartBinCheck(), FileBlockName::BinCheck_Len());
        uint32 nN=1;
        uint32 nBit;
        uiint nByte= sizeof(uiint);
        if(nByte==4) nBit=32;
        if(nByte==8) nBit=64;
        ofs.write((char*)&nN, 4);
        ofs.write((char*)&nBit, 4);
        ofs.write(FileBlockName::End(), FileBlockName::End_Len());
        oRes.WriteAlgebra_bin(ofs);
        oRes.WriteRes_bin(ofs);
    } else {
        ofs.open(sFileName.c_str(), ios::out);
        oRes.WriteAlgebra(ofs);
        oRes.WriteRes(ofs);
    }
    ofs.close();
}
void CFileWriterChunk::PrintResult_Start(const uiint& nStep, string filename, bool bBinary)
{
    stringstream ss;
    ss << nStep;
    string sFileName = filename + "." + ss.str();
    if(!mb_fstr) sFileName += ".rlt";
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::MWDebug,"Result Filename => ",sFileName);
    if(bBinary) {
        m_ofs.open(sFileName.c_str(), ios::out|ios::app|ios::binary);
    } else {
        m_ofs.open(sFileName.c_str(), ios::out|ios::app);
    }
    m_ofs << FileBlockName::StartResult() << endl;
}
void CFileWriterChunk::PrintResult(const uiint& width, string format, vector<void*>& param)
{
    iint *nVal;
    double* dVal;
    string sVal;
    uiint nLength = format.length();
    uiint ipara=0;
    for(uiint i=0; i < nLength; i++) {
        if(format[i] == '%') {
            ++i;
            switch( format[i] ) {
            case('d'):
                // width : 出力行の幅=一行の文字数
                nVal= (iint*)param[ipara];
                m_ofs << setw(width) << right << dec << *nVal << " ";
                ipara++;
                break;
            case('f'):
                dVal= (double*)param[ipara];
                m_ofs << setw(width) << right << fixed << *dVal << " ";
                ipara++;
                break;
            case('e'):
                dVal= (double*)param[ipara];
                m_ofs << setw(width) << right << scientific << *dVal << " ";
                ipara++;
                break;
            case('s'):
                sVal= (char*)param[ipara];
                m_ofs << setw(width) << right << sVal << " ";
                ipara++;
                break;
            default:
                break;
            }
        }
    };
    m_ofs << endl;
}
void CFileWriterChunk::PrintResult_End()
{
    m_ofs << FileBlockName::EndResult() << endl;
    m_ofs.close();
}
