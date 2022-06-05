% See Notes 2022-05-30-2100.
clear;
numFigs = 0;
setprngstates(0);
%
sizeF = 1;
sizeX = 50;
sizeU0 = 20;
numElemPerCol = 0;
numAddtlElemPerRow = 4;
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

matU0 = randn(sizeX,sizeU0);
matV = utorthdrop(randn(sizeX,sizeU0));
if (0)
	matU3 = zeros(sizeX,3);
	matU3(1:3:end,1) = 1.0;
	matU3(2:3:end,2) = 1.0;
	matU3(3:3:end,3) = 1.0;
	matU5 = zeros(sizeX,5);
	matU5(1:5:end,1) = 1.0;
	matU5(2:5:end,2) = 1.0;
	matU5(3:5:end,3) = 1.0;
	matU5(4:5:end,4) = 1.0;
	matU5(5:5:end,5) = 1.0;
	matV = [ matV, matU3, matU5 ];
endif
sizeV = size(matV,2);
%assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
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
if (0)
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
matResCollective = zeros(sizeF,sizeV);
for m=1:sizeF
	%sizeL = 50;
	sizeL = 7;
	assert( sizeL < sizeV );
	[ foo, orderedList ] = sort( matR(m,:) );
	elemUsed = orderedList(1);
	while (numel(elemUsed)<sizeL)
		sparsePrecon0__internal;
		
		elemUsed = [ elemUsed, newElemUsed ];
		
		%elemUsed = [ elemUsed, orderedList(1), orderedList(2) ]
		
		if (0)
		newElem = orderedList(foo<0.8);
		if ( 0 == numel(newElem) )
			break;
		endif
		elemUsed = [ elemUsed, newElem ];
		if ( numel(elemUsed) > sizeL )
			elemUsed = elemUsed(1:sizeL);
			break;
		endif
		endif
	endwhile
	%
	elemUsed
	usedMsk = logical(zeros(1,sizeX));
	usedMsk(elemUsed) = true;
	matWUsed = matW(m,:); % Actually just a row vector.
	matVUsed = matV(elemUsed,:);
	coeffs = matWUsed*(matVUsed')/(matVUsed*(matVUsed'));
	%
	matJCollectiveEst(m,elemUsed) = coeffs;
	%
	%
	matResCollective(m,:) = matW(m,:) - coeffs*matVUsed;
endfor
%matJCollectiveEst
%matJ

if (1)
	numFigs++; figure(numFigs);
	plot( matJ, 'o-', 'linewidth', 3, matJCollectiveEst, 'x-' );
	grid on;
	numFigs++; figure(numFigs);
	plot( matResCollective, 'o-' );
	grid on;
endif
