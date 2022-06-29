function [ varargout ] = Struct2mat(S)
%
fields = fieldnames(S);
nF = length(fields);
for iF=1:nF
    name = fields{iF};
    assignin('caller',name,S.(name));
    varargout{iF}=name;
end
end

