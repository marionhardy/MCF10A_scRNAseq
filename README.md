# MCF10A_scRNAseq
Single cell RNA seq data of MCF10A testing routine conditions. A companion dataset to MCF10A_bulkRNAseq.

------------------------------------ Where to find data and description ---------------------------------------------

All the data is stored in \\albecklab.mcb.ucdavis.edu\data\Notebooks\RNA Sequencing Project 2022-23\

These figures were produced using the code in .\Marion_Analysis\code\Oscar_scenic_v2_plotting.R
The data input is the csvb table in .\Oscar_Analysis\120723_scRNAseq_and_bulk\scRNAseq_OUTS.zip\scRNAseq_OUTS\SCENIC\albeck_pyscenic_analysis_v2_condition_specific_regulons_zscores.csv
The regulon enrichment data and z-scores were produced by Oscar Davalos through pySCENIC. This is his second iteration.
The figures were produced by me because Oscar's figures were too stringently filtered and didn't have clustering.

An extensive overview of the experimental set up can be found in .\Cell culture help chart_11.25.2022_MP.xlsx


--------------------------------------- Brief summary ----------------------------------------

MCF10A cells cultured for 16 hours in:

akt_inhibitor_ipasertip		IM spiked with an AKT inhibitor, ipasertib
ampk_activator_mk8722		IM spiked with an AMPK inhibitor, MK8722
base_imaging_medium    		base imaging medium (IM) (insulin, ct, hc, glucose, pyruvate)
erk_inhibitor_pd_0325901   	IM spiked with PD0325901, an ERK inhibitor                  
growth_medium			complete growth medium DMEMF12 (FBS, EGF, glutamine, hc, ct, insulin, glucose, pyruvate)                         
il6                       	IM with IL6
im_ct                      	IM without cholera toxin (CT)
im_egf                          IM with EGF     
im_glucose                	IM without glucose
im_glutamine                    IM without glutamine   
im_hc                  		IM without hydrocortisone (HC)
im_insulin                	IM without insulin
ldh_inhibitor_galloflavin 	IM spiked with a LDH inhibitor, galloflavin
m_torc1_inhibitor_rapamycin     IM spiked with a mTORC1 inhibitor, rapamycin  
mpc_inhibitor_uk5099            IM spiked with a MPC inhibitor, UK5099     
oligomycin 			IM spiked with an ETC inhibitor, oligomycin
                
------------------------------------------- NB --------------------------------------------

This data set comes from the scRNAseq data.
We also have live-cell images and bulk RNAseq available, from the same experiment.





























