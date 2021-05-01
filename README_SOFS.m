% README gives instructions on preparing data to run WDM.
%
% Wave staff data are simulated by computing x,y tilt from east-west and north-south buoy velocities.
% First the orbital velocities are inferred from the buoy velocities via
% the transfer function calculated in Orbital_vs_Buoy_Velocity_SOFS.m and
% stored in buoy_orbital_velocity_transfer_fcn.mat. Then tilt is obtained
% from the orbital velocity (increased by exp(k*edepth)) using linear
% theory; edepth is the effective action depth of the drag and inertia
% forces. These two steps are performed in the function:
% BuoyVelocityToSlope.m called by WDM.
%
% Direction estimates are good -- depending on the ratio of phase differences.
% Wavenumber estimates -- depending on the phase differences -- are good when the displacements 
% reflect the wave tilts, but at low frequencies (swell),gentle tilts, the current (non-wave) and wind 
% field dominate the horizontal displacements. Swell and Wind-sea are fine for both f
% & k directional spectra. Very low swell seems to be fine for the f-directional
% spectra, but is meaningless in the wavenumber spectra.
%
%
% Create 2 folders "fspect" and "kspect" to accept the spectra from
% "WAVEPLOTS" and "WAVENUM_LOG_SOFS" respectively. They are called from WDM.m 
% run by run. Create folder "wdms" to accept output of time series of
% wavenumbers, directions and amplitudes at each frequency from WDM.
%
% Run WDM. Assign run nos -- these will be used throughout the analyses.
% The programs are described in their headers.
% There are 2 streams and folders must be created in advance:
% Frequency stream:
%_____________________
% Program                    Write to Folder            What's written to folder.
% FITSECH_SOFS              beta                       mean directions and beta vs f.
% f_spreading_adjust_SOFS   fspectd                    corrected spectra, FCC(f,phi).
% kvsf_SOFS                 ---------                  plots the dispersion relation.
% f_polar_SOFS              ---------                  plots polar spectra.
%
% Wavenumber stream:
%_____________________
% Program                    Write to Folder            What's written to folder.
% FITSECHK_SOFS             betak                      mean directions and betak vs k.
% k_spreading_adjust_SOFS   kspectd                    corrected spectra, FCC(k,phi).
%                            NOTE:   f_spreading_adjust_SOFS must have been run before.
% k_polar_SOFS              ---------                  plots polar spectra.