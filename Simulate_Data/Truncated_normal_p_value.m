function p_val= Truncated_normal_p_value(max_stat,z_lo,varargin)
%Given an observation x from a truncated normal distibution (where the
%underlying normal has mean zero and variance one), 
%and z_lo and z_up are the lower and
%upper truncation points, respectively), this file calculates a one-sided
%p-value for the null hypothesis that the mean of x (again, in the
%underlying normal distribution) is zero, against the alternative that it
%is larger than zero.
%If z_up is not provided, it is assumed to be infinity

if( isempty(varargin) == 0)
    z_up = varargin{1};
else
    z_up = Inf;
end

if( ~ (z_lo <= max_stat && max_stat <= z_up) )
    warning('max_stat is not between z_lo and z_up');
end

%Adjust upper and lower values to improve numerical performance
V_lo_trim=max(z_lo,max_stat-10);
V_up_trim=min(z_up,max_stat+10);

f=@(y) exp(-.5*(y.^2-V_lo_trim.^2));
CDF_numerator=integral(f,V_lo_trim,max_stat,'RelTol',10^-6);
CDF_denominator=integral(f,V_lo_trim,V_up_trim,'RelTol',10^-6);

if CDF_denominator == Inf || CDF_denominator == -Inf || CDF_numerator == Inf || CDF_numerator == -Inf
    warning('Infinite integral in p-value calc')
end

p_val=1-CDF_numerator/CDF_denominator;

end

