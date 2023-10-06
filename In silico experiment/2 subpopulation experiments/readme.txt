/**
 * Functions and scripts need to generate the 2 subpopulation experiments and plots (Figures 3,4,5,6,7)
 *
 *  'main_PhenoPop_sto_2_pop.m': main script for 2 subpopulation experiments.
 *      Inside the script, we have 3 main sections:
 *          1. Initialization 
 *          2. Performing Bootstrapping
 *          3. Performing Point estimation
 *
 *      Expected outcome: 
 *          The results are saved as 'CI'+seed number+'(point_estimate).mat'.
 *
 *  'Plot_CI_EP_LC.m': Script for generating Illustrative example (Figure 3).
 *  
 *  'plot_CI_width.m': Script for generating the confidence interval width comparison plot (Figure 5).
 *
 *  'plot_pe_error.m': Script for generating the point estimation accuracy plots:
 *      Inside the script, we have 3 main sections:
 *          1. Record point estimation error from all 30 experiments.
 *          2. Plot the point estimation error of all modeling parameters (Figure 6).
 *          3. Plot the point estimation error of all important parameters (Figure 4).
 *  
 *  'Plot_joint_confidence_region.m': Script for generate the Joint confidence region of two parameters (Figure 7).
 *
 */
