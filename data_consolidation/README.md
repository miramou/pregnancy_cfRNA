# File description:

## Main files:
* merge_clean_data.R
  * **Important Note:**
    * To run this script, you must change the variable **mainDir** to reflect the path to the main directory as described in the main instructions
  * Pulls RT-qPCR data from csv files located in the folder **raw_data**
  * Linearly transform RT-qPCR threshold cycle values to estimated transcript counts and normalizes data to inputted plasma and qPCR volumes
  * Creates dataframe where samples from Denmark have been averaged using 3-week centered moving average
  * Saves data and gene lists of interests

## Auxillary files:
* calc_centered_moving_avg.R
  * Helper functions to calculate 3-week centered moving average
* calc_counts_per_mL.R
  * Helper functions to calculate estimated transcript counts per mL plasma
* combine_ga_dCD.R
  * Helper functions to merge metadata with existing dataframe
