# vIRt-Moore

### Simplified network model of ventral intermediate reticular zone neurons for whisking behavior

Last modified 2026-02-01

***

#### Instructions for sample simulations

- For the first simulation, run **virt_moore.m** to use default parameters
- For subsequent simulations, you can call the graphical user interface with **virt_moore_gui** and import parameter files that were generated


#### Instructions for reproducing each figure in Lu et al, 2026 Neuron ("Theory of amplitude control by frequency detuning in the rodent whisking oscillator circuit")

MATLAB Version used was **R2025a** and/or **R2025b** (results were identical for at least Figures 3 and 5). All required scripts for each figure are organized 
in it's respective subdirectory.

- **Figure 1**:
  - Run **run_rlc_simulations.m** in a desired output directory. This would create a subdirectory called **RLC_Sim_Test** that is timestamp labelled
- **Figure 2**:
  - Run **virt_analyze_sniff_whisk.m** in a directory containing a data subdirectory called **data_sniff_whisk**. This would create a subdirectory called **output_sniff_whisk** that contains all the outputs
- **Figure 3**:
  - Run **virt_moore.m** on its own using default parameters. This should generate the figures similar to 3e and 3f. Then call **virt_moore_gui** and import the last parameter file that was generated. Change the parameter **Square Wave Input -> pre-Botzinger Neurons -> Amplitude** to 0, then click **Run Simulation**. This should generate the figures 3c and 3d.
- **Figure 4a to 4m**:
  - Run **run_Lu_2026_Fig4ato4m_twister.m**
- **Figure 4n**:
  - Run **run_Lu_2026_Fig4n_threefry.m**
- **Figure 5**:
  - Run **run_Lu_2026_Fig5_threefry.m**
- **Figure 6 and Figure 7**:
  - Run **run_Lu_2026_Fig6and7_twister.m**