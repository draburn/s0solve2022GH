% Function...
 
function [ ary3Kappa, datOut ] = calcFDKappa( vecX, vecF, funchF, matV, matC, prm=[] )
	fevalCount = 0;
	%
	sizeX = size( vecX, 1 );
	assert( isrealarray(vecX,[sizeX,1]) );
	sizeF = size( vecF, 1 );
	assert( isrealarray(vecF,[sizeF,1]) );
	sizeV = size( matV, 2 );
	assert( isrealarray(matV,[sizeX,sizeV]) );
	assert( isrealarray(matC,[sizeV,sizeV]) );
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
	ary3Kappa = zeros( sizeX, sizeX, sizeF );
	% NOTE XXF NOT FXX!
	%
	%matJ = zeros(sizeF,sizeX);
	for n=1:sizeV
	if ( matC(n,n) )
		vecV = matV(:,n);
		vecXP = vecX; vecXP += epsX * vecV;
		vecXM = vecX; vecXM -= epsX * vecV;
		vecFP = funchF(vecXP); fevalCount++;
		vecFM = funchF(vecXM); fevalCount++;
		%matJ += ( vecFP - vecFM ) * (vecV') / (2.0*epsX);
		vecK = (vecFP+vecFM-(2.0*vecF))/(epsX^2);
		for m=1:sizeF
			ary3Kappa(:,:,m) += vecK(m)*(vecV*(vecV'));
		endfor
	endif
	endfor
	%
	for n1=1:sizeV
	for n2=1:n1-1
	if ( matC(n1,n2) || matC(n2,n1) )
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
		vecK = ( vecFPP + vecFMM - vecFPM - vecFMP ) / (4.0*(epsX^2));
		for m=1:sizeF
			ary3Kappa(:,:,m) += vecK(m)*(vecV1*(vecV2'));
			ary3Kappa(:,:,m) += vecK(m)*(vecV2*(vecV1'));
		endfor
	endif
	endfor
	endfor
	%
	if ( nargout >= 2 )
		datOut.fevalCount = fevalCount;
		%datOut.matJ = matJ;
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
%!	ary3Kappa_secret = randn(sizeX,sizeX,sizeF);
%!	for m=1:sizeF
%!		ary3Kappa_secret(:,:,m) = ( ary3Kappa_secret(:,:,m)' + ary3Kappa_secret(:,:,m) ) /2.0;
%!	endfor
%!	ary3Kappa_secret
%!	%
%!	funchF = @(x) funcVecQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret);
%!	vecX0 = randn(sizeX,1);
%!	vecF0 = funchF(vecX0);
%!	%
%!	sizeV = sizeX;
%!	matV = eye(sizeX,sizeV);
%!	matC = ones(sizeV,sizeV);
%!	ary3Kappa0 = calcFDKappa( vecX0, vecF0, funchF, matV, matC )
%!	assert( isrealarray(ary3Kappa0,[sizeX,sizeX,sizeF]) );
%!	assert( reldiff(ary3Kappa0,ary3Kappa_secret) < eps^0.3 );
%!	%
%!	sizeV = sizeX;
%!	matV = orth(randn(sizeX,sizeV));
%!	matC = ones(sizeV,sizeV);
%!	ary3Kappa0 = calcFDKappa( vecX0, vecF0, funchF, matV, matC )
%!	assert( isrealarray(ary3Kappa0,[sizeX,sizeX,sizeF]) );
%!	assert( reldiff(ary3Kappa0,ary3Kappa_secret) < eps^0.3 );
%!	%
%!	sizeV = 2;
%!	matV = eye(sizeX,sizeV);
%!	matC = eye(sizeV,sizeV);
%!	ary3Kappa0 = calcFDKappa( vecX0, vecF0, funchF, matV, matC )
%!	assert( isrealarray(ary3Kappa0,[sizeX,sizeX,sizeF]) );
%!	%
%!	sizeV = 2;
%!	matV = orth(randn(sizeX,sizeV));
%!	matC = eye(sizeV,sizeV);
%!	ary3Kappa0 = calcFDKappa( vecX0, vecF0, funchF, matV, matC )
%!	assert( isrealarray(ary3Kappa0,[sizeX,sizeX,sizeF]) );
