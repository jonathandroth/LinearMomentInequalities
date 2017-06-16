function p_val= Truncated_normal_p_value(max_stat,V_lo,varargin)
%Given an observation x from a truncated normal distibution (where the
%underlying normal has mean zero and variance one), 
%and V_lo and V_up are the lower and
%upper truncation points, respectively), this file calculates a one-sided
%p-value for the null hypothesis that the mean of x (again, in the
%underlying normal distribution) is zero, against the alternative that it
%is larger than zero.
%If V_up is not provided, it is assumed to be infinity

if( isempty(varargin) == 0)
    V_up = varargin{1};
else
    V_up = Inf;
end

%Adjust upper and lower values to improve numerical performance
V_lo_trim=max(V_lo,max_stat-10);
V_up_trim=min(V_up,max_stat+10);

f=@(y) exp(-.5*(y.^2-V_lo_trim.^2));
CDF_numerator=integral(f,V_lo_trim,max_stat,'RelTol',10^-6);
CDF_denominator=integral(f,V_lo_trim,V_up_trim,'RelTol',10^-6);
p_val=1-CDF_numerator/CDF_denominator;

end

