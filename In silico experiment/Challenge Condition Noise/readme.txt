/**
 * Functions and scripts need to generate the Challenge Condition noise experiments and plots (Figures 10,11,12,13)
 *
 *  'main_CC_Noise_boot.m': main script for Challenge Condition noise experiments bootstrapping results.
 *      Inside the script, we have 2 main sections:
 *          1. Initialization 
 *          2. Performing Bootstrapping
 *
 *      Expected outcome: 
 *          The results are saved as 'sup_Noise_quant'+seed number+'(var_fixed).mat' in Result folder.
 *
 *  'main_CC_Noise_point_estimate.m': main script for Challenge Condition noise experiments point estimation results.
 *      Inside the script, we have 2 main sections:
 *          1. Initialization 
 *          2. Performing point estimation
 *          
 *
 *      Expected outcome: 
 *          The results are saved as 'sup_Noise_quant'+seed number+'(point_estimation).mat' in Result folder.
 *
 *  'Plot_CI_no.m': Script for generating Illustrative example (Figure 10).
 *  
 *  'plot_CI_width.m': Script for generating the confidence interval width comparison plot (Figure 12,13).
 *
 *  'plot_quant_NO.m': Script for generating the point estimation accuracy plots (Figure 11).
 *
 */