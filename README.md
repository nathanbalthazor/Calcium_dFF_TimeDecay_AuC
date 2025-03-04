# Calcium_dFF_TimeDecay_AuC
A MATLAB script to analyze the differential fluorescent response of Calcium Indicator GCaMP to stimulation

This script was originally developed by Nathan Balthazor (ntbalthazor@gmail.com), a Research Assistant for Dr. Beth Habecker in the Chemical Physiology and Biochemistry department of Oregon Health Sciences University (OHSU). The Habecker lab focuses on autonomic imbalances in the heart, specifically exessive reinnervation of sympathetic efferents to the heart following cardiac injury such as cardiac infarction. 

This script was developed to analyze the data gathered by Dr. Arianna Scalco in her investigation of cardiac injury produced by hypertension. Dr. Scalco conducted calcium fluorescence imaging on mice expressing GCaMP6 in the neurons of the sympathetic stellate ganglia, half WT for ChAT and the other half with ChAT knocked-out in these cells, discouraging the formation of intraganglion cholinergic collaterals. Two experimental groups were created within these genotype populations, one control group and the other administered AngII via back pump, leading to hypertension. 

Stellate ganglia were surgically isolated from each mouse and placed in cell culture, where the culture bath was injected with Nicotine followed by KCl as a positive control. The resulting fluoresence was recorded by a Nikon camera using the Nikon NIS Elements software to export the data to MATLAB. In MATLAB, a github toolbox titled "Detect" by user nsdesai performed CNMF analysis on the video files to return a dF/F trace of each identifed neuron. 

This MATLAB data file of individual neuron dFF traces is able to be loaded into the code in this repository. By downloading TimeDecay_AuC_Algorithm and it's dedicated function file TimeDecay_AuC, large neuron x dFF .m files can be evaulated for time decay and area under the curve kinetics displayed by the dFF responses. The results are automatically plotted and saved, as well as numerical values exported to an excel file directly in the indicated directory. 

-----All questions for use or suggested improvements can be emailed to ntbalthazor@gmail.com----------
