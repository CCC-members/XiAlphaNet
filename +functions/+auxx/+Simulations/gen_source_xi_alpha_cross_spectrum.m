
options=init_simulation();


[~,theta0]     = gen_hggm02(Nseed,options);
[sigma,theta]  = spectral_precicion_tensor(theta0,options.nu0(2),...
options.tstudent_a(2),options.tstudent_b(2),options.tstudent_dof(2),...
options.deltaf,options.Fmin,options.Fmax);
process_energy = sigma_to_band_energy((options.Ntpoints)*sigma,...
    [options.band(2,1),options.band(2,2)],options.deltaf);