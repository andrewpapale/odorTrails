function Y = zScoreExperimenter(X,experimenter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nE = length(unique(experimenter));
Y = nan(size(X));
for iE=1:nE
    k = experimenter==iE;
    Y(k)=nanzscore(X(k));
end

end
