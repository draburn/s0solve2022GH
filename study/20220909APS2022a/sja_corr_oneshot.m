function [ matJ, datOut ] = sja_corr_oneshot( matV, matW, prm=[] )
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	maxNumNZEPerRow = mygetfield( prm, "maxNumNZEPerRow", floor(sizeK-sqrt(sizeK)) );
	tol = mygetfield( prm, "tol", 0.001 );
	[ foo, vecPriorityList ] = sort( diag(1.0./(eps+sumsq(matW,2))) * ((matW * (matV')).^ 2) * diag(1.0./(eps+sumsq(matV,2))), 2, "descend" );
	vecNZEList = vecPriorityList(:,1:maxNumNZEPerRow);
	matJ = zeros( sizeF, sizeX );
	for m = 1 : sizeF
		matJ(m,vecNZEList(m,:)) = (matV(vecNZEList(m,:),:)') \ (matW(m,:)');
	endfor
	vecRes = sum( (matJ*matV - matW).^2, 2 );
	vecNorm = sum( matW.^2, 2 );
	%semilogy( vecRes./vecNorm, 'o-' ); grid on;
	msg( __FILE__, __LINE__, sprintf( "Capture percent: %0.3f.", sum( vecRes <= (tol^2) * vecNorm )*100.0/sizeF ) );
	datOut = [];
return;
endfunction

% The first argument to "sort()" above is "matChi",
% which is calculated compactly above but can also be calculated as...
%for m=1:sizeF
%for n=1:sizeX
%	sumWV = 0.0;
%	sumWW = 0.0;
%	sumVV = 0.0;
%	for k=1:sizeK
%		w = matW(m,k);
%		v = matV(n,k);
%		sumWV += w*v;
%		sumVV += v*v;
%		sumWW += w*w;
%	endfor
%	matChi(m,n) =  (sumWV^2)/((eps+sumVV)*(eps+sumWW));
%endfor
%endfor
% See ztest_corr.m for more information.
