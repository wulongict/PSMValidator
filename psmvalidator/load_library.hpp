//
// Created by wulong on 4/5/17.
//
#include "../SpectraST/SpectraSTCreateParams.hpp"
#include "../SpectraST/SpectraSTLib.hpp"
#include "../SpectraST/SpectraSTLog.hpp"

// Verbose and quiet option flags. Don't want to pass them everywhere, so use global variables
bool g_verbose;
bool g_quiet;

// global pointer to a log object for keeping a log file
SpectraSTLog *g_log;

void generate_pepnovo_mgf_from_library()
{
    g_log = new SpectraSTLog("spectrast.log");
    string filename = "/data/wulong/data/NIST/test_new.splib";
    vector<string> fileNames;
    fileNames.push_back(filename);

    SpectraSTSearchParams searchParams;
    searchParams.indexCacheAll=true;

    SpectraSTLib *lib = new SpectraSTLib(filename, &searchParams,true);

    cout << "[Info] finished library " << lib << " contains entries" << endl;
    vector<SpectraSTLibEntry *> hits;
    lib->retrieve(hits,0,1202);

    cout << "[Info] found entries: " << hits.size() << endl;
    cout << "[Info] finished library " << lib << " contains entries"  << endl;
    SpectraSTPeakList * pkl = hits[0]->getPeakList();
    int num_peaks = pkl->getNumPeaks();  // you can not get the peptide sequence from pkl
    cout << "[Info] number of peaks in this entry " << num_peaks << endl;

    SpectraSTPeptideLibIndex *pepLibindex = lib->getPeptideLibIndexPtr();
    int pepnum = pepLibindex->getEntryCount();
    cout << "[Info] found peptide: " << pepnum << endl;
    pepLibindex->printStats(cout);


    delete lib;


}


void load_library()
{
    g_log = new SpectraSTLog("spectrast.log");
    string filename = "/data/wulong/data/NIST/test_new.splib";
    vector<string> fileNames;
    fileNames.push_back(filename);

    SpectraSTSearchParams searchParams;
    searchParams.indexCacheAll=true;

    SpectraSTLib *lib = new SpectraSTLib(filename, &searchParams,true);

    cout << "[Info] finished library " << lib << " contains entries" << endl;
    vector<SpectraSTLibEntry *> hits;
    lib->retrieve(hits,0,1202);

    cout << "[Info] found entries: " << hits.size() << endl;
    cout << "[Info] finished library " << lib << " contains entries"  << endl;
    SpectraSTPeakList * pkl = hits[0]->getPeakList();
    int num_peaks = pkl->getNumPeaks();  // you can not get the peptide sequence from pkl
    cout << "[Info] number of peaks in this entry " << num_peaks << endl;

    SpectraSTPeptideLibIndex *pepLibindex = lib->getPeptideLibIndexPtr();
    int pepnum = pepLibindex->getEntryCount();
    cout << "[Info] found peptide: " << pepnum << endl;
    pepLibindex->printStats(cout);


    delete lib;


}