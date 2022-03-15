% Function...
 
function [ matVTKV, datOut ] = calcFDVTKV( vecX, vecF, funchF, matV, prm=[] )
	assert( nargout <= 2 );
	assert( nargin <= 5 );
	fevalCount = 0;
	%
	sizeX = size( vecX, 1 );
	assert( isrealarray(vecX,[sizeX,1]) );
	sizeF = size( vecF, 1 );
	assert( isrealarray(vecF,[sizeF,1]) );
	sizeV = size( matV, 2 );
	assert( isrealarray(matV,[sizeX,sizeV]) );
	%
	doTests = true;
	if ( doTests )
		assert( reldiff(matV'*matV,eye(sizeV)) < eps^0.5 );
	endif
	%
	epsX = mygetfield( prm, "epsX", eps^0.25 );
	assert( isrealscalar(epsX) );
	assert( 0.0 < epsX );
	%
	%
	%
	% We could let the user pass in matJV, but,
	%  we'd want to be sure to avoid hitting same epsilon.
	% We'll consider only second-order finite differencing.
	matJV = zeros(sizeF,sizeV);
	ary3VTKappaV = zeros(sizeF,sizeV,sizeV);
	for n=1:sizeV
		vecV = matV(:,n);
		vecXP = vecX; vecXP += epsX * vecV;
		vecXM = vecX; vecXM -= epsX * vecV;
		vecFP = funchF(vecXP); fevalCount++;
		vecFM = funchF(vecXM); fevalCount++;
		matJV(:,n) = (vecFP-vecFM)/(2.0*epsX);
		ary3VTKappaV(:,n,n) = (vecFP+vecFM-(2.0*vecF))/(epsX^2);
	endfor
	%
	includeOffDiagonalTerms = mygetfield( prm, "includeOffDiagonalTerms", false );
	if ( includeOffDiagonalTerms )
		for n1=1:sizeV
		for n2=1:n1-1
			vecV1 = matV(:,n1);
			vecV2 = matV(:,n2);
			vecXPP = vecX; vecXPP += epsX*vecV1; vecXPP += epsX*vecV2;
			vecXPM = vecX; vecXPM += epsX*vecV1; vecXPM -= epsX*vecV2;
			vecXMP = vecX; vecXMP -= epsX*vecV1; vecXMP += epsX*vecV2;
			vecXMM = vecX; vecXMM -= epsX*vecV1; vecXMM -= epsX*vecV2;
			vecFPP = funchF(vecXPP); fevalCount++;
			vecFPM = funchF(vecXPM); fevalCount++;
			vecFMP = funchF(vecXMP); fevalCount++;
			vecFMM = funchF(vecXMM); fevalCount++;
			ary3VTKappaV(:,n1,n2) = ( vecFPP + vecFMM - vecFPM - vecFMP ) / (4.0*(epsX^2));
		endfor
		endfor
		for n1=1:sizeV
		for n2=n1+1:sizeV
			ary3VTKappaV(:,n1,n2) = ary3VTKappaV(:,n2,n1);
		endfor
		endfor
	endif
	%
	matVTKV = zeros(sizeV,sizeV);
	for n=1:sizeV
		matVTKV(n,:) += vecF'*ary3VTKappaV(:,:,n);
	endfor
	matVTKV = ( matVTKV' + matVTKV ) / 2.0;
	%
	if ( nargout >= 2 )
		datOut.fevalCount = fevalCount;
		datOut.matJV = matJV;
		datOut.matVTHV = matJV'*matJV + matVTKV;
		datOut.ary3VTKappaV = ary3VTKappaV;
		if ( mygetfield( prm, "calcSBKappaA", false ) );
			% Calculate standard basis kappa approximate.
			datOut.ary3SBKappaA = zeros(sizeF,sizeX,sizeX);
			for nf=1:sizeF
				datOut.ary3SBKappaA(nf,:,:) = matV * reshape( ary3VTKappaV(nf,:,:), [sizeV,sizeV] ) * (matV');
			endfor
		else
			datOut.ary3SBKappaA = [];
		endif
	endif
return;
end


%!test
%!	setprngstates(0);
%!	sizeX = 3;
%!	sizeF = 3;
%!	vecX_secret = randn(sizeX,1);
%!	vecF_secret = randn(sizeF,1);
%!	matJ_secret = randn(sizeF,sizeX);
%!	%
%!	ary3Kappa_secret = randn(sizeF,sizeX,sizeX);
%!	for nx1=1:sizeX
%!	for nx2=1:nx1-1
%!		ary3Kappa_secret(:,nx1,nx2) = (ary3Kappa_secret(:,nx1,nx2)+ary3Kappa_secret(:,nx2,nx1))/2.0;
%!	endfor
%!	endfor
%!	for nx1=1:sizeX
%!	for nx2=nx1+1:sizeX
%!		ary3Kappa_secret(:,nx1,nx2) = ary3Kappa_secret(:,nx2,nx1);
%!	endfor
%!	endfor
%!	%
%!	funchF = @(x) funcVecQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret);
%!	vecX0 = randn(sizeX,1);
%!	vecF0 = funchF(vecX0);
%!	%
%!	prm = [];
%!	prm.includeOffDiagonalTerms = true;
%!	prm.calcSBKappaA = true;
%!	[ matK, datOut ] = calcFDVTKV( vecX0, vecF0, funchF, eye(sizeX,sizeX), prm  );
%!	assert( isrealarray(matK,[sizeX,sizeX]) );
%!	assert( issymmetric(matK) );
%!	assert( isrealarray(datOut.matJV,[sizeF,sizeX]) );
%!	assert( isrealarray(datOut.ary3VTKappaV,[sizeF,sizeX,sizeX]) );
%!	%ary3Kappa_secret
%!	%datOut.ary3VTKappaV
%!	reldiff(ary3Kappa_secret,datOut.ary3VTKappaV)
%!	assert( reldiff(ary3Kappa_secret,datOut.ary3VTKappaV) < eps^0.33 );
%!	assert( reldiff(datOut.ary3VTKappaV,datOut.ary3SBKappaA) < eps^0.5 );
%!	%
%!	sizeV = ceil(sizeX/2);
%!	matV = orth(randn(sizeX,sizeV));
%!	prm = [];
%!	prm.includeOffDiagonalTerms = true;
%!	prm.calcSBKappaA = true;
%!	[ matVTKV, datOut ] = calcFDVTKV( vecX0, vecF0, funchF, matV, prm );
%!	assert( isrealarray(matVTKV,[sizeV,sizeV]) );
%!	assert( issymmetric(matVTKV) );
%!	assert( isrealarray(datOut.matJV,[sizeF,sizeV]) );
%!	assert( isrealarray(datOut.ary3VTKappaV,[sizeF,sizeV,sizeV]) );
%!	%
%!	for nv1=1:sizeV
%!	for nv2=1:sizeV
%!		vecFoo = zeros(sizeF,1);
%!		for nx1=1:sizeX
%!		for nx2=1:sizeX
%!			vecFoo += ary3Kappa_secret(:,nx1,nx2)*matV(nx1,nv1)*matV(nx2,nv2);
%!		endfor
%!		endfor
%!		assert( reldiff(datOut.ary3VTKappaV(:,nv1,nv2),vecFoo) < eps^0.33 )
%!	endfor
%!	endfor


%!test
%!	setprngstates(0);
%!	sizeX = 3;
%!	sizeF = 3;
%!	vecX_secret = randn(sizeX,1);
%!	vecF_secret = randn(sizeF,1);
%!	matJ_secret = randn(sizeF,sizeX);
%!	%
%!	ary3Kappa_secret = zeros(sizeF,sizeX,sizeX);
%!	for n=1:sizeX
%!		ary3Kappa_secret(:,n,n) = randn(sizeF,1);
%!	endfor
%!	funchF = @(x) funcVecQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret);
%!	vecX0 = randn(sizeX,1);
%!	vecF0 = funchF(vecX0);
%!	%
%!	sizeV = ceil(sizeX/2);
%!	matV = eye(sizeX,sizeV);
%!	[ matVTKV, datOut ] = calcFDVTKV( vecX0, vecF0, funchF, matV );
%!	assert( isrealarray(matVTKV,[sizeV,sizeV]) );
%!	assert( issymmetric(matVTKV) );
%!	assert( isrealarray(datOut.matJV,[sizeF,sizeV]) );
%!	assert( isrealarray(datOut.ary3VTKappaV,[sizeF,sizeV,sizeV]) );
%!	%
%!	for nv1=1:sizeV
%!	for nv2=1:sizeV
%!		vecFoo = zeros(sizeF,1);
%!		for nx1=1:sizeX
%!		for nx2=1:sizeX
%!			vecFoo += ary3Kappa_secret(:,nx1,nx2)*matV(nx1,nv1)*matV(nx2,nv2);
%!		endfor
%!		endfor
%!		assert( reldiff(datOut.ary3VTKappaV(:,nv1,nv2),vecFoo,eps) < eps^0.33 )
%!	endfor
%!	endfor
