function [] = qqplot_with_pvalue(x)

qqplot(x);
[~,pValue] = swtest(x);
pValue = round(pValue,3);

annotation('textbox', [.2 .5 .3 .3], 'String', strcat('Shapiro Wilk p-value:', num2str(pValue)), 'FitBoxToText','on' );

end