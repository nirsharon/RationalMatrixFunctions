% script name: run_all_figures   
% 
% All figures of Section 3, not included matrix functions
%
% NS, August 19

tic

% spline figures
figure_SPLINE_varing_bound
error_denominator_script

% other scalar test cases
figure_oscilatory_function
figure_cusp
figure_abs_val

% matrix function figures
figure_filter
figure_filter_vec_timing
figure_SPD_projection
figure_spd_timing

toc()