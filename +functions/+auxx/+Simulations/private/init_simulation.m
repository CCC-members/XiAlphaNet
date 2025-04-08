function options=init_simulation()

%% gen_hgggm options
options.config       = 2;
options.var          = 2;
options.connections  = [1 2; 2 1];
options.extensions   = [ceil(Nseed/3); ceil(Nseed/3); Nseed - 2*ceil(Nseed/3)];

% for spectrum waveform
Fs                   = 200;
Fmin                 = 0;
Fmax                 = Fs/2;
Ntpoints             = 512;
deltaf               = Fs/Ntpoints;
F                    = Fmin:deltaf:Fmax;
Nfreqs               = length(F);
options.Fs           = Fs;
options.Fmin         = Fmin;
options.Fmax         = Fmax;
options.F            = F;
options.Ntpoints     = Ntpoints;
options.deltaf       = deltaf;
options.Nsegments    = Nsegments;
options.Nfreqs       = Nfreqs;
% for xi and alpha process
options.nu0          = [0,10];
options.tstudent_a   = [900,600];
options.tstudent_b   = [5,9];
options.tstudent_dof = [3.2,60];
options.band         = [0,3;8,12];
% options.band         = [0,3;0,100];

% for xi process
options.vertices     = vertices;
options.faces        = faces;
options.elec_pos     = elec_pos;
%
options.Nsource      = Nsource;
options.Nsensor      = Nsensor;
options.Nsubj        = Nsubj;
options.Nseed        = Nseed;
options.Nsim         = Nsim;


end



