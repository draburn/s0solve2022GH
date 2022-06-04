% See Notes 2022-05-30-2100.
clear;
numFigs = 0;
setprngstates(0);
%
sizeF = 1;
sizeX = 5;
sizeV0 = 4;
numElemPerCol = 0;
numAddtlElemPerRow = 2;
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
matV = utorthdrop(randn(sizeX,sizeV0));
sizeV = size(matV,2);
assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
matW = matJ*matV;
%
matWVAvg = matW*(matV') / sizeV;
vecWSqAvg = sum( matW.^2, 2 ) / sizeV;
vecVSqAvg = sum( matV.^2, 2 ) / sizeV;
%
matJIndivEst = zeros(sizeF,sizeX);
for n=1:sizeX
	matJIndivEst(:,n) = matWVAvg(:,n) / vecVSqAvg(n);
endfor
%matJIndivEst
%
matWVAvgSq = matWVAvg.^2;
matR = zeros(sizeF,sizeX);
for m=1:sizeF
for n=1:sizeX
	matR(m,n) = 1.0 - ( matWVAvgSq(m,n) / ( eps + vecWSqAvg(m)*vecVSqAvg(n) ) );
endfor
endfor
%matR
%
if (1)
	m = 1;
	numFigs++; figure( numFigs );
	plot( matJ(m,:), 'o-' );
	grid on;
	numFigs++; figure( numFigs );
	plot( matJIndivEst(m,:), 'o-' );
	grid on;
	numFigs++; figure( numFigs );
	plot( 1.0./matR(m,:), 'o-' );
	grid on;
	vecR = matR(m,:);
	vecA = 1.0./vecR;
	vecA = vecA .* ( vecR <= sort(vecR)(sizeV-1) );
	numFigs++; figure( numFigs );
	plot( vecA, 'o-' );
	grid on;
	%return;
	matJ
	matV
	matW
	matJIndivEst
	matR
endif
%
%
%
% Use quick method: just pick top several.
useQuickMethod = false;
if ( useQuickMethod )
	for m=1:sizeF
		sizeL = sizeV-1;
		[ foo, orderedList ] = sort( matR(m,:) );
		elemSelector = orderedList(1:sizeL)
	endfor
	error( "Quick method not implemented; already looks flawed on first trial." );
endif
%
%
%
matJCollectiveEst = zeros(sizeF,sizeX);
for m=1:sizeF
	sizeL = 2;
	[ foo, orderedList ] = sort( matR(m,:) );
	elemUsed = orderedList(1)
	while (numel(elemUsed)<sizeL)
		sparsePrecon0__internal;
		elemUsed = [ elemUsed, newElemUsed ];
	endwhile
	%
	usedMsk = logical(zeros(1,sizeX));
	usedMsk(elemUsed) = true
	matWUsed = matW(m,:) % Actually just a row vector.
	matVUsed = matV(usedMsk,:)
	coeffs = matWUsed*(matVUsed')/(matVUsed*(matVUsed'))
	%
	matJCollectiveEst(m,usedMsk) = coeffs;
endfor
matJCollectiveEst
matJ
