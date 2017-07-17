function p_val= Truncated_normal_p_value_by_simulation(max_stat,z_lo,varargin)
%Given an observation x from a truncated normal distibution (where the
%underlying normal has mean zero and variance one), 
%and z_lo and z_up are the lower and
%upper truncation points, respectively), this file calculates a one-sided
%p-value for the null hypothesis that the mean of x (again, in the
%underlying normal distribution) is zero, against the alternative that it
%is larger than zero.
%If z_up is not provided, it is assumed to be infinity

%Computation is done using 20,000 draws from the truncated normal using the
%algorithm of:
% Z. I. Botev (2015), "The Normal Law Under Linear Restrictions:
%  Simulation and Estimation via Minimax Tilting", submitted to JRSS(B)
if( isempty(varargin) == 0)
    z_up = varargin{1};
else
    z_up = Inf;
end

if( ~ (z_lo <= max_stat && max_stat <= z_up) )
    warning('max_stat is not between z_lo and z_up');
end

%For whatever reason, mvrandn doesn't work with a one-dimensional variable,
% so we draw from a 2-dimensional normal without any correlation and then
% stack the 2 rows
rng(0);
draws = mvrandn([z_lo;z_lo],[z_up;z_up],eye(2),10000);
draws = draws(:);

p_val = mean( max_stat <= draws); %pvalue is fraction of draws greater than maxtat


end



