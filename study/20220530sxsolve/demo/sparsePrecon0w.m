% See Notes 2022-05-30-2100.
clear;
numFigs = 0;
setprngstates(0);
%
sizeF = 1;
sizeX = 20;
numElemPerCol = 0;
numAddtlElemPerRow = 3;
c0 = 0.0;
csx = 0.0;
csf = 0.0;
%
matJ0 = zeros(sizeF,sizeX);
for n=1:sizeX
for t=1:numElemPerCol
	m = ceil( sqrt(eps) + (sizeF-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = randn();
endfor
endfor
for m=1:sizeF
for t=1:numAddtlElemPerRow
	n = ceil( sqrt(eps) + (sizeX-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = randn();
endfor
endfor
matJ0 += c0 * randn(sizeF,sizeX);
assert( isrealarray(matJ0,[sizeF,sizeX]) );

matSX = diag(exp(csx*randn(sizeX,1)));
matSF = diag(exp(csf*randn(sizeF,1)));
matJ = matSF*matJ0/matSX;



%setprngstates(96551984); % Fail
%setprngstates(75868048); % Succ
%setprngstates();
sizeU0 = 15;

matU0 = randn(sizeX,sizeU0);
matV = utorthdrop(randn(sizeX,sizeU0));

sizeV = size(matV,2);
%assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
matW = matJ*matV;



matV2 = matV.^2;
if (1)
	matWeight = 1.0./(eps+1.0-matV2); % Weights.
else
	matWeight = ones(sizeX,sizeV);
endif
matJEst = zeros(sizeF,sizeX);
matCEst = zeros(sizeF,sizeX);
matChi = zeros(sizeF,sizeX);
matRho = zeros(sizeF,sizeV);
for m=1:sizeF
	%
	rvecW = matW(m,:);
	rvecCEst = zeros(1,sizeX);
	rvecJEst = zeros(1,sizeX);
	elemUsed = [];
	newElem = [];
	for l=1:5
		% This is not quite right; these shouldn't be frozen.
		rvecJEst(newElem) = rvecCEst(newElem);
		rvecRho = rvecW - rvecJEst*matV;
		%
		%rvecWVAvg = sum( rvecW.*(matV.*matWeight), 2 ); % Automatic broadcastig.
		%rvecV2Avg = sum( matV2.*matWeight, 2 );
		%rvecW2Avg = sum( (rvecW.^2).*matWeight, 2 ); % Automatic broadcasting.
		rvecRVAvg = sum( rvecRho.*(matV.*matWeight), 2 ); % Automatic broadcastig.
		rvecV2Avg = sum( matV2.*matWeight, 2 );
		rvecR2Avg = sum( (rvecRho.^2).*matWeight, 2 ); % Automatic broadcasting.
		rvecCEst = rvecRVAvg./rvecV2Avg;
		rvecChi = (rvecRVAvg.^2)./( eps + rvecR2Avg.*rvecV2Avg );
		%
		%%%[ foo, orderedList ] = sort( -rvecChi );
		%%%elemUsed = orderedList(1:l) % Already "used" elements will remain used.
		usedElemMsk = logical(zeros(1,sizeX));
		usedElemMsk(elemUsed) = true;
		selectorList = (1:sizeX)(~usedElemMsk);
		[ foo, orderedList ] = sort( -rvecChi(~usedElemMsk) );
		newElem = selectorList(orderedList(1));
		elemUsed = [ elemUsed, newElem ]
	endfor
	vecJEst(elemUsed) = rvecCEst(elemUsed);
	cvecRho = rvecW - rvecJEst*matV;
	%
	matJEst(m,:) = rvecJEst;
	matCEst(m,:) = rvecCEst;
	matChi(m,:) = rvecChi;
	matRho(m,:) = rvecRho;
endfor
%
m = 1;
numFigs++; figure(numFigs);
plot( matJ(m,:), 'o-', 'linewidth', 5, matCEst(m,:), 's-', 'linewidth', 1, matJEst(m,:), 'x-', 'linewidth', 2 );
grid on;
numFigs++; figure(numFigs);
plot( matChi(m,:), '^-', 'linewidth', 2 );
grid on;
numFigs++; figure(numFigs);
plot( matRho(m,:), '^-', 'linewidth', 2 );
grid on;




return;
