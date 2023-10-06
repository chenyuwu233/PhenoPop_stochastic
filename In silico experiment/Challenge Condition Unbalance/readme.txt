/**
 * Functions and scripts need to generate the Challenge Condition unbalanced proportion experiments and plots (Figures 14,15)
 *
 *  'main_CC_Unbalance.m': main script for Challenge Condition unbalanced proportion experiments bootstrapping results.
 *      Inside the script, we have 2 main sections:
 *          1. Initialization 
 *          2. Performing Bootstrapping
 *
 *      Expected outcome: 
 *          The results are saved as 'Close_GR'+seed number+'(var_fixed).mat' in Result folder.
 *
 *  'main_CC_Unbalance_estimate.m': main script for updating the point estimation results.
 *      Inside the script, we have 2 main sections:
 *          1. Initialization 
 *          2. Performing point estimation
 *          
 *
 *      Expected outcome: 
 *          The results are saved as 'Close_GR'+seed number+'(point_estimation).mat' in Result folder.
 *
 *  'Plot_CI_hp.m': Script for generating Illustrative example (Figure 14).
 *
 *  'plot_quant_hp.m': Script for generating the point estimation accuracy plots (Figure 15).
 *
 */