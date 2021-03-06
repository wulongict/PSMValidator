//
// Created by wulong on 11/4/18.
//

#ifndef MYTOOL_DATABASEMANAGER_H
#define MYTOOL_DATABASEMANAGER_H

#include <sqlite3.h>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <map>
using namespace std;
class CTable;

class CDBEntry{
    map<string, vector<string>> m_data;
    bool m_colname_fixed;
public:
    int size(){
        return m_data.begin()->second.size();
    }
    bool empty(){return m_data.empty() or m_data.begin()->second.empty();}
    CDBEntry(){m_colname_fixed = false;}
    CDBEntry(const vector<string>& headers);
    void setValue(string header, int record_id, string value){
        if(m_data.count(header)==0 or record_id >= m_data[header].size()){
            throw runtime_error("Error: header and record out of range");
        }
        m_data[header][record_id] = value;
    }
    void extendTable(string header){
        if(m_data.find(header)!=m_data.end()){
            // found new header already exist. do nothing.
        }else{
            int n = size();
            m_data[header]=vector<string>(n,"");
            // fill in blank strings
        }
    }
    // add data
    void add(const string& header, const string& value);
    string get(string header, int record_id){
        string result;
        if(m_data.find(header)==m_data.end()){
            cout << "Invalid header" << endl;
        }
        if(m_data[header].size()<=record_id) {
            cout << "Invalid column --" << header << endl;
            cout << "valid header are " << endl;
            for(auto &x: m_data){
                cout << x.first << endl;
            }
        }
        result=m_data[header].at(record_id);
        return result;
    }
    long getInt(string header, int record_id);
    double getFloat(string header, int record_id);

    void print();
};


// Here is how a callback function looks like
//typedef int (*sqlite3_callback)(
//        void*,    /* Data provided in the 4th argument of sqlite3_exec() */
//        int,      /* The number of columns in row */
//        char**,   /* An array of strings representing fields in the row */
//        char**    /* An array of strings representing column names */
//);

static int get_row_callback(void *data, int numCol, char **argv, char **azColName){
    // the columns of a gt table
    // ID = 0
    //FILEID = 0
    //MS2COUNTS = 0
    //PEPTIDE = UNKNOWN
    //SCORE = -1.0
    //SCAN = 11
    //CTERM = 0.0
    //NTERM = 0.0
    //MODIFICATION = UNMODIFIED
    //PRECURSOR = 462.151
    //CHARGE = 1
    //RT = 5.42103
    //PEPTIDEPROPHETPROB = 0.0
    //IPROPHETPROB = 0.0
    //RFSCORE = 0.0
    //ISDECOY = 0
    //SIGNIFICANCE = 0
    //PROTEIN =
    //CE =
    //ALTERPEPTIDE = NULL
    //NEIGHBOR =
    int i;
    vector<string> *pdata = static_cast<vector<string>*>(data);
    pdata->assign(numCol,"");
    for(i=0; i<numCol; i++){
//        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");

        (*pdata)[i] = argv[i]?argv[i]:"NULL";
    }
    //printf("\n");
    return 0;
}


static int get_row_callback_dbEntry(void *data, int numCol, char **argv, char **azColName){
    CDBEntry *pdata = static_cast<CDBEntry*>(data);

    for(int i=0; i<numCol; i++){
        pdata->add(azColName[i],argv[i]?argv[i]:"NULL");
    }
    return 0;
}

// vector<vector<string> > * data 
static int get_multiple_rows_callback(void *data, int numCol, char **argv, char **azColName){
    int i;
    vector<vector<string>> *pdata = (vector<vector<string>>*)data;

    pdata->emplace_back(vector<string>(numCol,""));
    int row_num = pdata->size()-1;
    for(i=0; i<numCol; i++){
        //printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
        (*pdata)[row_num][i] = argv[i]? argv[i]:"NULL";
    }
    //printf("\n");
    return 0;
}



class CDataBaseManager
{
    sqlite3 *m_db;

    void openDB(const string& databaseFileName);
public:
    CDataBaseManager(const string& databaseFileName);

    ~CDataBaseManager();

    void getRow(vector<string> &result, const string& idname, int idvalue, const string& tablename, bool verbose);
    void getRow_deprecated(vector<string> &result, const string& sql, bool verbose);
    void getRow(CDBEntry &result, const string& sql, bool verbose);
    void getMultipleRows(vector<vector<string>> &result, const string& sql, bool verbose);


    void execAsTransaction(const string& sql, bool verbose);

    bool tableExists(const string& tablename, bool verbosity);

    long getTotalRows(string tablename);

    void batchSQL(bool run, vector<string> &manySQLs, const string& sql, int batchsize, bool verbose);


};


void ConvertTable();


// This two functions are use to create the sql table from begining. However, we may want to update it gradually.
// SpecFiles table
// id, fileid, start, end
void batchInsertSpecfiles(CTable &specfiles, CDataBaseManager &cbm, int start, int end);

// groundtruth table
// id, fildid, ms2counts, peptide, score, scan, cterm, nterm, modification, precursor, charge, rt, peptide prophet prob, iprophet prob.
void batchInsertGroundTruth(CTable &groundtruth, CDataBaseManager &cbm, int start, int end);

#endif //MYTOOL_DATABASEMANAGER_H
