//
// Created by wulong on 11/25/17.
//
// Nov 30. fixed the

#include "../librarymsms/Util.h"
#include "boost/program_options.hpp"
#include "dependentlibs/msbasic/ConcretePSMFeatures.h"
#include "dependentlibs/msbasic/CDebugMode.h"
#include "dependentlibs/msbasic/CFlow.h"
#include "spdlog/spdlog.h"
#include "dependentlibs/msbasic/commonfun.h"

void generateConfigFile(const boost::program_options::options_description &patternLR, string filename);

boost::program_options::variables_map getParam(int argc, char *argv[]) {
    namespace po = boost::program_options;
    po::options_description patternLR("Command line only"), configs_infile("Config file");
    patternLR.add_options()
            ("config", po::value<string>(),
             "the config file")
            ("help,h", po::bool_switch()->default_value(false), "print help information");
    configs_infile.add_options()
           ("inputfile,i", po::value<string>(), "input file[mgf generated by SpectraST]")
            ("task", po::value<string>(),
             "[pepxml2feat|psmfeature|trainfragmodel|fragscore|cometsearch]:  "
             "psmfeature. extract psm feature. "
             "trainfragmodel. train model for fragmentation pattern. "
             "fragscore. get fragmentation score alone. "
             "cometsearch. cometsearch")
            // PSM feature workflow =========================
            ("cid", po::bool_switch()->default_value(false),
             "mass accuracy. cid for low mass accuracy, by default it's high mass accuracy, that is cid=false")
            ("annotation,a",
             po::value<string>()->default_value("searchresult"),
             "annotation method of input mgf file;")
            ("searchresult", po::value<string>()->default_value(""),
             "specify the comet search result if the annotation method is searchresult")

            ("debug", po::bool_switch()->default_value(false),
             "only work on user specified psm . turn on verbosity; only for psm feature.")
            ("debugOutputNullPep", po::bool_switch()->default_value(false),
             "output null peptide .")
            ("debugUseScanCharge", po::bool_switch()->default_value(false),
             "use scan+charge to look for ans .")
            ("debugRank", po::value<int>()->default_value(0), "rank of search result to use ")

            // Liblinear model: The peakpair model
            ("use_output", po::bool_switch()->default_value(false),
             "export training and testing file for liblinear executables.")
            ("min_foldchange,m", po::value<double>()->default_value(2),
             "minimal intensity fold change to select matched peak pairs as samples ")
            ("model", po::value<string>()->default_value(""), "filename of pre-trained model")
            ("peptidelen", po::bool_switch()->default_value(false), "build model for different peptide length")
            ("ghostpeak,g", po::bool_switch()->default_value(false),
             "use fragment ion with zero intensity in pattern learning")
            ("overwrite", po::bool_switch()->default_value(false), "redo all the steps and overwrite any file.")
            ("scoretype,s", po::value<string>()->default_value("counts"),
             "scoring method for combining scores of individual peaks in each MS2. [prob2|logprob|prob|weightedprob|counts]")

            // Comet Search parameters ======================
            ("database", po::value<string>()->default_value(""), "database for comet search")
            ("targetdb", po::value<string>()->default_value(""), "database for target-only comet search")
            ("decoydb", po::value<string>()->default_value(""),  "database for decoy-only comet search")
            ("outname", po::value<string>()->default_value(""),  "output name for comet search")

            ("cometparam", po::value<string>()->default_value("/data/wulong/data/proteometools/comet.params.new"),
             "comet parameter file")
            ("cometbinary", po::value<string>()->default_value("/tools/bin/comet.2016012.linux.exe"), "comet binary path")
            ("specfile", po::value<string>()->default_value(""), "input spectra file")
            ("decoytype,d", po::value<string>()->default_value("decoy"),
             "how to find decoy file? decoy type [decoy|secondhit]. decoy --> *decoyfrag.txt secondhit --> *secondhitfrag.txt")

            ("truthfile", po::value<string>()->default_value(""), "the ground truth of psm")
            ("useAlternativeProt", po::value<bool>()->default_value(true),
             "to fix protein annotation of a peptide, use alternative protein list to reannotate. Priority on Target instead of decoy. default: false ")
            ("truthmethod", po::value<string>()->default_value("pepset"),
             "use a list of peptide (pepset) or a list of scan-peptide pairs (scan+pep) as groundtruth")

            ("trainValidator", po::bool_switch()->default_value(false),
             "train Validator (true) or use existing model (false); Default: false.")
            ("validatorProb", po::bool_switch()->default_value(false),
             "train Validator with probability estimation (true) or binary output (false); Default: true.")
            ("validatorModel", po::value<string>()->default_value("psmvalidatorRFmodel.foresst"),
             "the name of model used for calculate probability of each PSM.")
            ("featurelistfile", po::value<string>()->default_value("/data/wulong/data/NIST/comet2016012/validator/featurelist.txt"),
             "the feature list file .")
            ("writeValidatorModel", po::value<string>()->default_value("psmvalidatorRFmodel.forest"),
             "the file name of validator model; used when trainValidator parameter is set as true; the file extension must be .forest")
            ("featureXtractionHitRank", po::value<int>()->default_value(0),
             "the hit ranking: 0,1,2,3,4 used for feature extraction. default: 0, the top ranking peptide is used.")

            // The Ranger: Random Forest for PSM score ============================================
            ("rangerbinary", po::value<string>()->default_value("/data/wulong/bitbucket/codejam/tools/cmake-build-release/Release/ranger"),
             "ranger binary path")
            ("trainingSampleSize", po::value<int>()->default_value(90000),
             "Number of samples used for training. Input feature file must contain more feature than needed.) ")
            ("mtry", po::value<int>()->default_value(3),
             "Number features for split, sqrt(total_feature) ")
            ("ntree", po::value<int>()->default_value(200),
             "number of trees/estimators in the forest ")
            ("maxDepthAllowed", po::value<int>()->default_value(6),
             "Maximum tree depth allowed in the forest ")
            ("thread", po::value<int>()->default_value(0),
             "Number of threads to be used. Zero means using #(cpu)-1 threads ")
            ("verbose,v", po::bool_switch()->default_value(false), "more detail output about progress on screen");

    patternLR.add(configs_infile);
    generateConfigFile(configs_infile,File::CFile(argv[0]).basename+"_default_param.conf");

    boost::program_options::variables_map vm;
    po::positional_options_description p;
    p.add("inputfile", -1);
    po::store(po::command_line_parser(argc, argv).options(patternLR).positional(p).run(), vm);

    if (vm.count("config")) {
        std::ifstream ifs{vm["config"].as<string>().c_str()};
        if (ifs) {
            po::store(parse_config_file(ifs, configs_infile), vm);
        }
    }
    po::notify(vm);

    if (vm.count("help") and vm.at("help").as<bool>()) {
        cout << "Build: " << __DATE__ << " " << __TIME__ << endl;
        cout << patternLR << endl;
        throw invalid_argument("Program Finish");
    }
    displayParam(vm);
    return vm;
}

void generateConfigFile(const boost::program_options::options_description &patternLR, string filename) {
    ofstream fout(filename.c_str(),ios::out);
    for(auto &eachOption: patternLR.options()){
        if("config"==eachOption->long_name() or "help" == eachOption->long_name()) {
            cout << "eachOption " << eachOption->long_name() << endl;
            continue;
        }

        string typeName;
        string defaultValue;

        boost::any value;
        eachOption->semantic()->apply_default(value);

        ostringstream  oss;
        if (typeid(float) == value.type()) {
            typeName = "float";
            oss <<  boost::any_cast<float>(value) ;
        }
        if (typeid(double) == value.type()) {
            typeName = "double";
            oss <<  boost::any_cast<double>(value) ;
        }
        if (typeid(string) == value.type()) {
            typeName = "string";
            oss <<  boost::any_cast<string>(value);
        }
        if (typeid(int) == value.type()) {
            typeName = "int";
            oss << boost::any_cast<int>(value) ;
        }
        if (typeid(bool) == value.type()) {
            typeName = "bool";
            oss <<  boolalpha << boost::any_cast<bool>(value) ;
        }
        defaultValue = oss.str();

        fout << "# " << eachOption->description() << endl
            << "# type: " << typeName << endl
             << eachOption->long_name() << " = "  << defaultValue << endl << endl;

    }
}

CFlow *CreateFlow(int argc, char *argv[]) {
    auto vm = getParam(argc, argv);
    string binaryPath = CPath(argv[0]).m_folder;
    cout << "Binary path is " << binaryPath << endl;
    string inputfile = vm.at("inputfile").as<string>();

    bool isSepPepLen = vm.at("peptidelen").as<bool>(); // not using
    bool use_ghost_peak = vm.at("ghostpeak").as<bool>();
    bool verbosity = vm.at("verbose").as<bool>();
    string scoretype = vm.at("scoretype").as<string>();
    double minIntensityFC = vm.at("min_foldchange").as<double>();
    string fragPatternFileName = vm.at("model").as<string>();

    bool use_output = vm.at("use_output").as<bool>();
    bool debug = vm.at("debug").as<bool>();
    bool overwrite = vm.at("overwrite").as<bool>();
    bool isCID = vm.at("cid").as<bool>();
    string taskname = vm.at("task").as<string>();
    string searchresult = vm.at("searchresult").as<string>();
    string annotationMethod = vm.at("annotation").as<string>();
    string featurelistfile = vm.at("featurelistfile").as<string>();

    // Comet Search
    string database = vm.at("database").as<string>();
    string outname = vm.at("outname").as<string>();
    string cometparam = vm.at("cometparam").as<string>();
    string cometbinary = vm.at("cometbinary").as<string>();
    string targetdb = vm.at("targetdb").as<string>();
    string decoydb = vm.at("decoydb").as<string>();

    string specfile = vm.at("specfile").as<string>();
    bool isTrainingRFModel = vm.at("trainValidator").as<bool>();
    bool isProbOutput = vm.at("validatorProb").as<bool>();
    string rfModelOutput = vm.at("writeValidatorModel").as<string>();
    string rfModelRead = vm.at("validatorModel").as<string>();
    int trainingSampleSize = vm.at("trainingSampleSize").as<int>();

    string truthFile = vm.at("truthfile").as<string>();
    string truthmethod = vm.at("truthmethod").as<string>();

    string rangerbinary = vm.at("rangerbinary").as<string>();
    int threadNum = vm.at("thread").as<int>();
    int mtry = vm.at("mtry").as<int>();
    int ntree = vm.at("ntree").as<int>();
    int maxDepth = vm.at("maxDepthAllowed").as<int>();
    int hitrank = vm.at("featureXtractionHitRank").as<int>();
    bool useAternativeProt = vm.at("useAlternativeProt").as<bool>();

    if (fragPatternFileName.empty()) {
        fragPatternFileName = isTrainingRFModel? getfragmodelfile(rfModelOutput, minIntensityFC, isTrainingRFModel): getfragmodelfile(rfModelRead, minIntensityFC, isTrainingRFModel);
    }

    spdlog::get("A")->info("Creating workflow.. ");
    CDebugMode *pDebug = CDebugMode::callDebug();
    pDebug->setDebug(debug);
    pDebug->outputNullPep = vm["debugOutputNullPep"].as<bool>();
    pDebug->rank = vm["debugRank"].as<int>();
    pDebug->useScanCharge = vm["debugUseScanCharge"].as<bool>();

    CFlow *ptr;
    if (taskname == "psmfeature") {
        ptr = new FeatureWorkFlow(inputfile, use_ghost_peak, minIntensityFC,
                                  fragPatternFileName, scoretype, annotationMethod, debug,
                                  "", featurelistfile, isCID, binaryPath);
    } else if (taskname == "trainfragmodel") {
        ptr = new FragModel(inputfile, verbosity, minIntensityFC, overwrite, fragPatternFileName, use_output,
                            use_ghost_peak, isCID, binaryPath);
    } else if (taskname == "fragscore") {
        ptr = new FragPatternScoreFlow(inputfile, verbosity, minIntensityFC, overwrite, fragPatternFileName, scoretype,
                                       searchresult, isCID, binaryPath);
    } else if (taskname == "cometsearch") {
        ptr = new CometSearch(cometbinary, database, cometparam, outname, inputfile);
    } else if (taskname == "cometsearchtd") {
        ptr = new CometSearchTD(cometbinary, targetdb, decoydb, cometparam, inputfile);
    } else if (taskname == "all") {
        ptr = new FlowAll(scoretype, use_ghost_peak, use_output, overwrite, minIntensityFC, verbosity,
                          fragPatternFileName, inputfile, cometbinary, targetdb, decoydb, cometparam,
                          isTrainingRFModel, isProbOutput, rfModelOutput, rfModelRead, mtry, ntree,
                          featurelistfile, isCID, trainingSampleSize, maxDepth, rangerbinary, binaryPath);
        // check parameter
        if (not isTrainingRFModel) {
            if (not File::isExist(rfModelRead)) {
                spdlog::get("A")->error("Random forest model file {} does not exist!", rfModelRead);
                throw invalid_argument(
                        "\nPlease try to add random forest model with option: --validatorModel XXX.forest\n");
            }
        }
    } else if (taskname == "ranger") { // testing ranger here
        // Input: CSV file with all the features.
        // e.g.
        // sepal_length,sepal_width,petal_length,petal_width,class
        // 5.1,3.5,1.4,0.2,1
        // 4.9,3.0,1.4,0.2,1
        // 4.7,3.2,1.3,0.2,1

        // input file should be specified as the feature.txt
        string feature_tsv_file;

        ptr = new RangerWraper(feature_tsv_file, isTrainingRFModel, rfModelRead, true, mtry, ntree, maxDepth, rangerbinary);
        // check parameter
        if (not isTrainingRFModel) {
            if (not File::isExist(rfModelRead)) {
                spdlog::get("A")->error("Random forest model file {} does not exist!", rfModelRead);
                throw invalid_argument(
                        "\nPlease try to add random forest model with option: --validatorModel XXX.forest\n");
            }
        }
    } else if (taskname == "tsv2rangerformat") {
        // Testing: to be deleted!
        string posFeature = "/data/wulong/data/proteometools/c57miss5decoy.pep.xmlfeatures.feat",
                negFeature = "/data/wulong/data/proteometools/test_new.mgffeatures.feat";
        string combinedPNfeature = inputfile + "_testing_PNcombined.feat";
        ptr = new RangerFormat(posFeature, negFeature, combinedPNfeature, trainingSampleSize);
    } else if (taskname == "validator") {
        shared_ptr<CTruth> truth = CreateTruth(truthFile, truthmethod);
        string pepxmlfilename = inputfile;
        // todo: support other format, like
        //  PepXML of : comet search result
        //  Sptxt: library spectra
        ptr = new ValidatePSM(pepxmlfilename, rfModelRead, scoretype, minIntensityFC, use_ghost_peak, threadNum,
                              truth, mtry,
                              ntree, featurelistfile, maxDepth, useAternativeProt, rangerbinary, binaryPath);
    } else if (taskname == "pepxml2feat") {

        string pepxmlfilename = inputfile;
        //
//        ptr = new ExtractFeatures(pepxmlfilename, rfModelRead, scoretype, minIntensityFC, use_ghost_peak, threadNum,
//                                  featurelistfile, useAternativeProt, binaryPath);
        ptr = new ExtractFeaturesFromPepXML(pepxmlfilename, rfModelRead, scoretype, minIntensityFC, use_ghost_peak,
                                            threadNum,
                                            featurelistfile, useAternativeProt, binaryPath, hitrank);
    }else if (taskname == "plot") {
        // in this case, the inputfile is pepxml file
        shared_ptr<CTruth> truth = CreateTruth(truthFile, truthmethod);
        cout << "[resultAnalysis] Attention: no basename " << endl;
        ptr = new resultAnalysis(inputfile, truth, "");
    } else {
        cout << "[Error] Incorrect taskname: \"" << taskname << "\"" << endl;
        ptr = nullptr;
    }
    return ptr;
}

void displayTitle() {
    spdlog::get("A")->info("\n"
                           "-------------------------------------------------\n"
                           "#            PSM validator                      #\n"
                           "#       Build Date: {} {}        #\n"
                           "#     Developed in Henry Lam's Group @ HKUST    #\n"
                           "-------------------------------------------------", __DATE__, __TIME__);
}

int main(int argc, char *argv[]) {
    initlog("psmvalidator.log", "A");
    displayTitle();
    spdlog::get("A")->info("CMD: {}", argToStr(argc, argv));

    try {
        SimpleTimer st("PSMValidator");
        CFlow *cf = CreateFlow(argc, argv);
        if (cf == nullptr) {
            spdlog::get("A")->error("Fail to create workflow! Program will exit!");
        } else {
            cf->run();
            delete cf;
            CDebugMode::releasePtr();
            spdlog::get("A")->info("Program finished successfully!");
        }
    }
    catch (const exception &ex) {
        cerr << "program exit with exception: " << ex.what() << endl;
    } catch (...) {
        cerr << "Unknown exceptions caused crash. Program terminated!" << endl;
    }
    spdlog::drop_all();
    return 0;
}

