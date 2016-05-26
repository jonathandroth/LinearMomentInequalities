function p_val= Truncated_normal_p_value(max_stat,V_lo)
%Given an observation x from a truncated normal distibution (where the
%underlying normal has variance one), and V_lo and V_up are the lower and
%upper truncation points, respectively), this file calculates a one-sided
%p-value for the null hypothesis that the mean of x (again, in the
%underlying normal distribution) is zero, against the alternative that it
%is larger than zero.

f=@(y) exp(-.5*(y.^2-V_lo.^2));
CDF_numerator=integral(f,V_lo,max_stat,'RelTol',10^-6);
CDF_denominator=integral(f,V_lo,max(10^5,max_stat+100),'RelTol',10^-6);
p_val=1-CDF_numerator/CDF_denominator;

end

