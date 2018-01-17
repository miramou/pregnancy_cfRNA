# pregnancy_cfRNA
Repository to reproduce analyses from manuscript titled "A noninvasive blood test for fetal development predicts gestational age and preterm delivery"

# To run the code

## Install the dependencies
### R packages
* plyr
* dplyr
* ggplot2
* reshape2
* cowplot
* corrplot
* caret
* randomForest
* edgeR
* limma
* locfit
* splines
* KernSmooth
* statmod
* MASS
* exactRankTests
* pROC
* effsize

### Python packages
* seaborn
* numpy
* pandas
* matplotlib

## Set up directory structure
The code expects the following directory structure:
* Main directory contains the following files and subdirectories:
  * Files:
    * load_data.R (Located in **common_scripts** folder)
    * plot_theme.R (Located in **common_scripts** folder)
    * counts_gene_lists.RData (Located in **common_scripts** folder)
  * Sub-directories titled as followed:
    * raw_data
      * Contains 9 csv files
    * data_consolidation
      * Contains 4 R scripts
    * fig_1
      * Contains 1 R script
    * fig_2
      * Contains 4 sub-sub-directories
    * fig_3
      * Contains 2 sub-sub-directories
    * supplement
      * Contains 5 R scripts

### Download raw data and relevant code
* Data can be found in the folder raw_data.
* All relevant code for each figure can be found in the folder with the corresponding name

### To run code
Please note that you must change the variable **mainDir** in any script to reflect the main directory (See section titled "*Set up directory structure*" 
#### From scratch:
* Download the raw data files and process them using the script **merge_clean_data.R** prior to running any other script.
* Build ML models using the script **build_models.R** in the folder **fig_2/c_d**
#### Using processed data:
* Download counts_gene_lists.RData located in the folder **common_scripts** prior to running any other script.
* Download saved model objects from the folders **fig_2/c_d/models** and **fig_2/e/models**
