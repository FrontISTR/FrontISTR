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
#include <unistd.h>
#include "HEC_MW3.h"
#include "Partitioner.h"
#include <getopt.h>

using namespace pmw;
using namespace std;

void usage()
{
    cout << "usage: " << "part" << " -i FILE -o FILE -n NUMBER" << endl;
    exit(0);
}

int main(int argc, char** argv)
{
    string inputPrefix, outputPrefix;
    uiint nparts = 0;

    //handle command line options
    struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"input", required_argument, NULL, 'i'},
        {"output", required_argument, NULL, 'o'},
        {"nparts", required_argument, NULL, 'n'},
        {0,0,0,0}
    };

    int result, option_index;
    while((result=getopt_long(argc, argv, "i:o:n:h", long_options, &option_index))!=-1) {
        switch(result) {
        case 'h':
            usage();
            break;
        case 'i':
            inputPrefix = optarg;
            break;
        case 'o':
            outputPrefix = optarg;
            break;
        case 'n':
            nparts = atoi(optarg);
            break;
        case '?':
            //usage();
            break;
        }
    }
    if(inputPrefix.empty()) usage();
    if(outputPrefix.empty()) usage();
    if(nparts <= 0) usage();

    //setup HECMW
    CMW *pMW = CMW::Instance();
    pMW->Initialize(argc, argv);
    uiint nDisplay= pMW->getDisplayDevice();
    uiint nInfo=pMW->getInfoMode(), nErr=pMW->getErrorMode(), nWarn=pMW->getWarnMode();
    pMW->LoggerMode(nInfo);
    pMW->LoggerDevice(nInfo, nDisplay);
    pMW->LoggerDevice(nErr,  nDisplay);
    pMW->LoggerDevice(nWarn, nDisplay);
    pMW->Banner_Display();

    // read mesh file
    pMW->FileRead(inputPrefix, ASCII);
    pMW->SelectAssembleModel(0);//lebel 0

    // setup neighbor data
    pMW->LoggerInfo(nInfo, "%s", "HEC_MW3 setupNeighbors");
    pMW->setupNeighbors();

    // partitioner
    CPartitioner p(pMW, nparts);
    pMW->LoggerInfo(nInfo, "%s", "part prepare");
    p.prepare();
    pMW->LoggerInfo(nInfo, "%s", "part partition");
    p.partition();
    pMW->LoggerInfo(nInfo, "%s", "part post process");
    p.post();
    pMW->LoggerInfo(nInfo, "%s", "part write");
    p.write(outputPrefix);

    p.printReport();;

    pMW->Finalize();
    return 0;
}
