function [dx,minMSE,nSelect] = dxdt(x,dT,window,postSmoothing)

% dx = dxdt(x,varargin)
% window = 1; % seconds
% postSmoothing = 0.5; % seconds --- 0 means don't
%
% Based on Janabi-Sharifi/Hayward/Chen, Discrete-time adaptive windowing
% for velocity estimation, IEEE Transactions on Control Systems Technology,
% (2000) 8(6):1003-1009.
% But modified extensively.  Basic algorithm is to allow windows from 3
% steps to nW = window/DT steps.  For each window, let dx = x(i+nW) - x(i).
% Select the window with the smallest MSE = sum_k=1..nW (x(i+k) - linear-fit(x(i+k) given slope from dx)).^2.
%
% postSmoothing does a convolution of normalized ones(nPS/DT)
%
% ADR 2003
% CORRECTED ADR 6 August 2012 - it was returning the negative direction

% window = 1; % seconds
% postSmoothing = 0.1; % seconds --- 0 means don't
display = 0;
%process_varargin(varargin);

% x = ctsd(removeNaNs(x));
% xD = x.data();
% dT = x.dt();
xD = x(~isnan(x(:)));

%dT = 1/50;

nW = min(ceil(window/dT),length(xD));
nX = length(xD);

MSE = zeros(nX, nW);
b = zeros(nX,nW);

MSE(:,1:2) = Inf;
nanvector = nan(nW,1);

for iN = 3:nW
	if display, fprintf(2,'.'); end
	b(:,iN) = ([nanvector(1:iN); xD(1:(end-iN))] - xD)/iN;
	for iK = 1:iN
		q = ([nanvector(1:iK); xD(1:(end-iK))] - xD + b(:,iN) * iK);
		MSE(:,iN) = MSE(:,iN) + q.*q;		
	end
	MSE(:,iN) = MSE(:,iN)/iN;	
end
if display, fprintf(2, '!'); end

[minMSE0, nSelect0] = min(MSE,[],2);
dx = nan .* ones(size(xD));
for iX = 1:nX
	dx(iX) = -b(iX,nSelect0(iX)) / dT;  % CORRECTED ADR 6 August 2012 - it was returning the negative direction
end

if postSmoothing
	nS = ceil(postSmoothing/dT);
	dx = conv2(dx,ones(nS)/nS,'same');
end

dx0 = nan(size(x));
dx0(~isnan(x(:)))=dx;

minMSE = nan(size(x));
nSelect = nan(size(x));
minMSE(~isnan(x(:)))=minMSE0;
nSelect(~isnan(x(:)))=nSelect0;

dx = dx0;
	
%dx = tsd(x.range(),dx);
