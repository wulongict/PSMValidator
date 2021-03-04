//
// Created by wulong on 3/13/21.
//

#include "gtest/gtest.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <vector>
#include <memory>
#include "CFlow.h"
#include "CDebugMode.h"

using namespace std;

TEST(ROC_,AUC_){
    initlog("x.log", "A");
    std::vector<double> positiveScores, negativeScores;
    positiveScores.push_back(0.1);
    negativeScores.push_back(0.2);
    positiveScores.push_back(0.2);
    negativeScores.push_back(0.3);
    cout << "Hi" << endl;
    CROCPlot roc( negativeScores,positiveScores),rocSwitch( positiveScores,negativeScores);
    cout << "hello" << endl;
    cout << rocSwitch.getAUC() << endl;
    cout << roc.getAUC() << endl;
    cout << "Hello" << endl;

}