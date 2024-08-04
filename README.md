# Matched-filtering and deboosting scripts
This repository contains the set of useful python scripts for the overdensity analysis on the SCUBA-2 850$\mu$m data, which is used for the paper "RAGERS I: Evidence of Anisotropic Distribution of
Submillimeter Galaxies in the 4C23.56 Field". However, these python scripts are easily adapted for other single-dish data.

Content of the main repository:
- crop_map.py: Crop the calibrated map into an appropriate size, in order to remove the noisy data on the outskirts.
- MF.py: Apply a matched filter (optimized for point sources) on the calibrated science map and do the source extraction through a simple clean-like algorithm.
- sim_proto.py: Generate 10,000 mock catalogs with a hard-coded overdensity. Use the blank field mock maps to estimate the overdensity is required.
- new_comp.py: Estimate the completeness factor and deboosting value from mock catalogs. 
- fidelity.py: Calculate the average number of spurious sources as a function of SNR, flux, or noise.
- statistic.py: Show how the deboosting, completeness, and fidelity factors vary as a function of SNR.
- offset.py: Estimate the position offsets as a function of SNR.
- deboost_dist.py: Deboost the catalogue following the Bayesian statistic for 10,000 times, which is used to estimate the protocluster number counts.
- nc_plot.py: Display the protocluster number counts and compare with the Geach+2017 and Stach+2018. 

Content in blank:
- s850_generate.py: Generate 10,000 mock catalogs based on the blank field number counts (Geach et al. 2017) with the same noise distribution as the calibrated data, which is used to estimate the overdensity of the protocluster field.

This repository is still under construction. We plan to publish on pypi and a jupyter notebook to demonstrate how our codes work soon.

Note: we realize that using "pure noise map" to obtain contamination rate will lead to a severe underestimation. It is because the presence of astronomical sources will also cause more spurious sources in the map. We have another treatment to overcome this issue. However, to have a fair comparison with Geach+17/Simpson+19, the current method is still suggested.

Please feel free to try our codes on other single-dish data of any protocluster field. We are curious to see how the SMG distribution of other protocluster field compare with the result in our paper. If you have any issue, please contact us at dazzh@space.dtu.dk
