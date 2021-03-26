
# PSMValidator
Validate PSM with a universal score.
## 1. Installation
### 1.1 Install with binary package
You could also simply try the compiled binary package. Download from [here](https://github.com/wulongict/PSMValidator/releases/tag/v1.0.0)  
Those are the binary files in the bin folder:  
```tree
.
|-- bin
|   |-- comet
|   |-- predict
|   |-- psmvalidator
|   |-- ranger
|   `-- train
```

### 1.2 Compile from source code
#### 1.2.1 Prerequisites
The following packages are required to compile from source code  

|package|version|provided|function|
|:-:|:-:|:-:|:-:|
|spdlog|x|YES|logger|
|comet|x|YES|  generate ground truth PSMs|
|ranger|x|YES|  random forest|
|liblinear|2.21|YES| Logistic Regression |
|rapidxml|1.13|YES|  xml parser|
|SpectraST|5.0|YES|spectra annotation|
|Boost|1.65|YES|File system; program option parser|
|gnuplot-iostream|x|YES|generate plots|
|sqlite3|3.0|NO||
|gnuplot|x|NO|needed by gnuplot-iostream|
|gsl|1.16|NO|SpectraST|

The following script installs the prerequisites packages.
```bash
# part one, install GSL 1.16
GSL_VERSION=$(gsl-config --version)
if [ "$GSL_VERSION" == "1.16" ]; then
  echo "GSL version 1.16 already installed!"
else
  echo GSL version is $GSL_VERSION
  wget https://mirror-hk.koddos.net/gnu/gsl/gsl-1.16.tar.gz
  tar xvf gsl-1.16.tar.gz
  # shellcheck disable=SC2164
  (
    cd gsl-1.16 || exit
    ./configure
    make all -j 10
    sudo make install
  )
fi

# part two install other packages.: sqlite3, gnuplot
sudo apt-get install sqlite3 libsqlite3-dev gnuplot gnuplot-qt
```


#### 1.2.2 Compilation and installation from source code (Ubuntu)
First, install cmake, build-essential and gnuplot with following commands
```bash
sudo apt-get install cmake build-essential gnuplot
```
Next, compile and install the psmvalidator with following commands.
```bash
unzip PSMValidator_v1.0.0.zip
cd </path/to/the/source/code>  # change accordingly
cmake .
cmake --build build
cmake --install . --prefix .
```

## 2 Basic usage
### 2.1 Tree structure after installed

📦build  
┣ 📂bin  
┃ ┣ 📜comet  
┃ ┣ 📜predict  
┃ ┣ 📜psmvalidator  
┃ ┣ 📜ranger  
┃ ┗ 📜train  
┣ 📂param  
┃ ┗ 📜psmvalidator.conf  
┣ 📂scripts  
┃ ┗ 📜train_small_demo_mgf.bash  
┗ 📂train  
┃ ┣ 📜comet16low.param  
┃ ┣ 📜featurelist.txt  
┃ ┣ 📜features.txt  
┃ ┣ 📜human_hcd_selected_new_small.mgf  
┃ ┣ 📜uniprot-human-2020-12.fasta  
┃ ┣ 📜uniprot-human-2020-12TD_only_decoy.fasta  
┃ ┣ 📜uniprot_yeast_reviewed_6721_Nov152016.fasta  
┃ ┗ 📜uniprot_yeast_reviewed_6721_Nov152016TD_only_decoy.fasta  

### 2.2 Train a demo validator
**input**:    
- spectra (mgf):   
    - train/train_small_demo_mgf                   
- database (fasta):  
    - target: uniprot-human-2020-12.fasta  
    - decoy:  uniprot-human-2020-12TD_only_decoy.fasta   
    
**output**:   
- model: 
    - human_hcd_selected_new_small_mtry_8_ntree_900_features.forest
      
```bash
cd scripts
./train_small_demo_mgf.bash
```

### 2.3 Prediction
**input**: 
- Comet search result (PepXML): 
- Model: human_hcd_selected_new_small_mtry_8_ntree_900_features.forest  

**output**:
- results table with RF (random forest) score: 

To use the validator, users could directly use the trained model. Here is an example:
```bash
cd scripts
./validate_psm.bash
```


