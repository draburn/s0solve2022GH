function [ matJ, datOut ] = sja_corr( matV, matW, prm=[] )
	mydefs;
	sizeX = size(matV,1);
	sizeF = size(matW,1);
	sizeK = size(matV,2);
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	maxNumNZEPerRow = mygetfield( prm, "maxNumNZEPerRow", sizeK-1 );
	tol = mygetfield( prm, "tol", sqrt(eps) );
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isscalar(valdLev) );
		assert( isscalar(verbLev) );
		assert( isposintscalar(sizeX) );
		assert( isposintscalar(sizeF) );
		assert( isposintscalar(sizeK) );
		assert( isrealarray(matV,[sizeX,sizeK]) );
		assert( isrealarray(matW,[sizeF,sizeK]) );
		assert( isposintscalar(maxNumNZEPerRow) );
		assert( maxNumNZEPerRow <= sizeK );
		assert( 0.0 < tol );
		assert( tol < 1.0 );
	endif
	%
	%
	% Try one-shot.
	tryOneShot = false;
	if ( tryOneShot )
	[ foo, vecPriorityList ] = sort( diag(1.0./(eps+sumsq(matW,2))) * ((matW * (matV')).^ 2) * diag(1.0./(eps+sumsq(matV,2))), 2, "descend" );
	vecNZEList = vecPriorityList(:,1:maxNumNZEPerRow);
	matJ = zeros( sizeF, sizeX );
	for m = 1 : sizeF
		matJ(m,vecNZEList(m,:)) = (matV(vecNZEList(m,:),:)') \ (matW(m,:)');
	endfor
	vecRes = sum( (matJ*matV - matW).^2, 2 );
	vecNorm = sum( matW.^2, 2 );
	oneShotTol = mygetfield( prm, "oneShotTol", 0.01 );
	oneShotCaptureRows = sum( vecRes <= oneShotTol^2 * vecNorm );
	msg( __FILE__, __LINE__, sprintf( "One-shot capture percentage: %0.3f (%0.3f required).", ...
	  oneShotCaptureRows*100.0/sizeF, sqrt(sizeF)*100.0/sizeF ) );
	if ( oneShotCaptureRows < sqrt(sizeF) )
		matJ = [];
		datOut = [];
		return;
	endif
	endif
	%
	%
	% DRaburn 2022-09-22:
	% Realized that identifying an element to be zero is a possibility;
	% really, it's not a question of number of Non-Zero-Elements per row,
	% but, rather, the number of determinable elements per row.
	% Something to keep in mind for next version!
	%
	matA = matV';
	%
	matJ = zeros(sizeF,sizeX);
	matE = 1+sparse(sizeF,sizeX); % 1 = undetermined element; 0 = determined elem "NZE".
	vecProcess = 1+zeros(sizeF,1); % Process this row?
	matNZEList = zeros(sizeF,0); % List of NZEs for each row.
	matD = matW; % Residual, "difference".
	doMainLoop = true;
	for numNZEPerRow = 1 : maxNumNZEPerRow
		vecProcess = sumsq( matD, 2 ) > (tol^2)*sumsq( matW, 2 );
		if ( 0 == sum(vecProcess) )
			break;
		endif
		%[ foo, vecNextNZE ] = max( abs(matD * matA) .* matE, [], 2 );
		% See below for explanation.
		[ foo, vecNextNZE ] = max( diag(1.0./(eps+sumsq(matD,2))) * ((matD * (matV')).^ 2) * diag(1.0./(eps+sumsq(matV,2))), [], 2 );
		matNZEList = [ matNZEList, vecNextNZE ];
		for nf = (1 : sizeF)(vecProcess)
			matE(nf,matNZEList(nf,end)) = 0;
			matJ(nf,matNZEList(nf,:)) = matA(:,matNZEList(nf,:)) \ (matW(nf,:)');
		endfor
		matD = matW - matJ*matV;
	endfor
	%
	datOut = [];
	return;
return;
endfunction

% The first argument to "max()" above is "matChi",
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
