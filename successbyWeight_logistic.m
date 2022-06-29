% logistic regression
X = [imw',pwgh'];
Xvars = {'','mouse','perweight'}';
Y = sfw';
options.MaxIter = 10000;
[b0,~,stats] = glmfit(X,Y,'binomial','link','logit','options',options);
table(Xvars,stats.p, b0, stats.p<0.05/2)

% ANOVAs
[~,bw] = histc(pwgh,quantile(pwgh,3));
[p,tab,stats] = anovan(nvw,{imw,bw},'model','interaction','varnames',{'mouse','per weight'});

[~,bw] = histc(pwgh,quantile(pwgh,3));
[p,tab,stats] = anovan(dnidw,{imw,bw},'model','interaction','varnames',{'mouse','per weight'});

[~,bw] = histc(pwgh,quantile(pwgh,3));
[p,tab,stats] = anovan(idpw,{imw,bw},'model','interaction','varnames',{'mouse','per weight'});

