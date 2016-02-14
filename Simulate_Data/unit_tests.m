

%% Test the sim_offeringss_firmf function
global T;
T = 10000;

p_o = .5;
p_d = .25;

[t,f, g, in_t , in_tminus1] = sim_offerings_firmf(1, p_o, p_d);

mean( in_t) %Should be close to stationary values p_o / (p_o + p_d) = 2/3 
mean( in_t( in_tminus1 == 1) ) %Should be close to 1 - p_d = .75
mean( in_t( in_tminus1 == 0) ) %Should be close to p_o = .5


%% Test simulated J for all firms

mean(in_t) * G * F %this should be approximately equal to numproducts

mean( in_t(in_tminus1 == 1)) %this should be aprpoximately 1 - p_d
mean( in_t(in_tminus1 == 0 & f == 1)) %this should be aprpoximately p_o_f(1)

