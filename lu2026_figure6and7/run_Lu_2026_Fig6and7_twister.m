% run_Lu_2026_Fig6and7_twister.m

jm_automate_virt_sim_vary_ext_current('RngAlgorithm', 'twister'); 
jm_automate_rerun_virt_sim_from_params_intrinsic_vs_driven_freq('RngAlgorithm', 'twister'); 
jm_postprocess_virt_sim_ExtCurrent_analysis; 
jm_postprocess_virt_sim_amplitude_analysis_ExtCurrentRuns;