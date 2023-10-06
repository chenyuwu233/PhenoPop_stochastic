/**
 * Functions and scripts need to generate the small time observation window experiments and plots (Figures 8)
 *
 *  'main_initial_fluctuation.m': main script for small time observation window experiments.
 *      Inside the script, we have 3 main sections:
 *          1. Initialization 
 *          2. Performing Bootstrapping estimation
 *          3. Performing Point estimation
 *
 *      Expected outcome: 
 *          The results are saved as 'CI_init_time_',seed number,'.mat' in the Result folder.
 *
 *  'plot_result.m': Script for generating small time observation results (Figure 8).
 *      This script contains 2 plots corresponding to accuracy plot and precision plot.
 *
 */