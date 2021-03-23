//
// Created by wulong on 3/15/21.
//

#ifndef PSMVALIDATOR_SABCDEFGH_H
#define PSMVALIDATOR_SABCDEFGH_H


#include <iostream>
#include <thread>
#include <algorithm>
#include <chrono>
#include <cfloat>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>
#include <numeric>
#include <vector>
#include <iterator>

using namespace std;


class CDBEntry;



struct SPsmAnnotation {
    // new members
    long idx;
    int ms2idx;
    int fileid;

    int isDecoy;
    int significance;
    string protein;

    string peptideseq;// = "UNKNOWN";
    double score;// = -1;
    double cterm_mass;// = 0;
    double nterm_mass;// = 0;
    string modificationstr;// = "UNMODIFIED";
    double precursormass;// = -1;
    int charge;// = -1;
    int ms2_scan;
    string mzxml_filename;
    double retentiontimeinsec;// = 0;
    double pProb;// = 0;
    double iProb;// = 0;
    double rfscore;//=0
    double precursorNeutralMass() const    {
        double pmass = 1.007276;
        return precursormass * charge - charge *pmass;
    }
    string m_collision_energy; // collision energy
    int peaknum;
    string m_neighbors;

    SPsmAnnotation();

    friend
    ostream &operator<<(ostream & out, const SPsmAnnotation & gt);


    SPsmAnnotation(int scan, double precursor, int chg, double rt);

    void set(int scan, double precursor, int chg, double rt, int fileId, long total_idx, string specfile, int ms2counts);
    void setCollisionEnergy(string collisionEnergy);

    void initWithRow(CDBEntry &results);

    string createInsertSql() const;
    string createUpdateSql() const;

    void outputToFile(ofstream &fout) const    {
        fout << idx << "\t" << mzxml_filename << "\t" << fileid << "\t" << ms2idx << "\t" << peptideseq
             << "\t" << score << "\t" << ms2_scan << "\t" << cterm_mass << "\t" << nterm_mass
             << "\t" << modificationstr << "\t" << precursormass << "\t" << charge << "\t" << retentiontimeinsec
             << "\t" << pProb << "\t" << iProb << "\t" << isDecoy << "\t" << significance << "\t" << protein << endl;
    }

    void toOstringStreamNoId(ostringstream &oss) const;

    void initMembers();

    string getModifiedPeptide(bool verbose) const;
    bool isSig()const {return significance == 1;}

    struct IdxDistPair{
        struct idxDist{
            double dist;
            long idx;
        };
        vector<idxDist> m_data;

    };


    IdxDistPair  neighborStrToIdxDistPair( char delimitor_1 = ';', char delimitor_2='@') const;
};
void getmodificationfrompeptidestring(string peptidestr, SPsmAnnotation &gtinfo);

#endif //PSMVALIDATOR_SABCDEFGH_H
