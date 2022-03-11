% Function...
%  Calculates the "Kappa array" (quadratic term) of a function via second-order finite-differencing.
 
function [ matJ0, ary3Kappa0, datOut ] = calcFDJKappa( vecX0, funchF, prm=[] )
	fevalCount = 0;
	%
	sizeX = size( vecX0, 1 );
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	vecF0 = mygetfield( prm, "vecF0", [] );
	if ( isempty(vecF0) )
		vecF0 = funchF( vecX0 ); fevalCount++;
	endif
	sizeF = size( vecF0, 1 );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	vecEpsFD = mygetfield( prm, "vecEpsFD", [] );
	if ( isempty(vecEpsFD) )
		epsFDRel = mygetfield( prm, "epsFDRel", eps^0.25 );
		assert( isrealscalar(epsFDRel) );
		assert( 0.0 < epsFDRel );
		vecEpsFD = epsFDRel * ( 1.0 + abs(vecX0) );
	endif
	assert( isrealarray(vecEpsFD,[sizeX,1]) );
	assert( sum(double(0.0<=vecEpsFD)) == sizeX );
	%
	%
	%
	matJ0 = zeros(sizeF,sizeX);
	ary3Kappa0 = zeros(sizeF,sizeX,sizeX);
	for n=1:sizeX
		epsFD = vecEpsFD(n);
		vecXP = vecX0; vecXP(n) += epsFD;
		vecXM = vecX0; vecXM(n) -= epsFD;
		vecFP = funchF(vecXP); fevalCount++;
		vecFM = funchF(vecXM); fevalCount++;
		matJ0(:,n) = (vecFP-vecFM)/(2.0*epsFD);
		ary3Kappa0(:,n,n) = (vecFP+vecFM-(2.0*vecF0))/(epsFD^2);
	endfor
	for n1=1:sizeX
	for n2=1:n1-1
		epsFD1 = vecEpsFD(n1);
		epsFD2 = vecEpsFD(n2);
		vecXPP = vecX0; vecXPP(n1) += epsFD1; vecXPP(n2) += epsFD2;
		vecXPM = vecX0; vecXPM(n1) += epsFD1; vecXPM(n2) -= epsFD2;
		vecXMP = vecX0; vecXMP(n1) -= epsFD1; vecXMP(n2) += epsFD2;
		vecXMM = vecX0; vecXMM(n1) -= epsFD1; vecXMM(n2) -= epsFD2;
		vecFPP = funchF(vecXPP); fevalCount++;
		vecFPM = funchF(vecXPM); fevalCount++;
		vecFMP = funchF(vecXMP); fevalCount++;
		vecFMM = funchF(vecXMM); fevalCount++;
		ary3Kappa0(:,n1,n2) = ( vecFPP + vecFMM - vecFPM - vecFMP ) / (4.0*epsFD1*epsFD2);
	endfor
	endfor
	for n1=1:sizeX
	for n2=n1+1:sizeX
		ary3Kappa0(:,n1,n2) = ary3Kappa0(:,n2,n1);
	endfor
	endfor
	%
	if ( nargout >= 2 )
		datOut.fevalCount = fevalCount;
		datOut.vecEpsFD = vecEpsFD;
		datOut.vecF0 = vecF0;
		datOut.vecG0 = matJ0'*vecF0;
		datOut.matH0 = matJ0'*matJ0;
		for nf=1:sizeF
			datOut.matH0 += reshape( vecF0(nf)'*ary3Kappa0(nf,:,:), [sizeX,sizeX] );
		end
	endif
return;
end


%!test
%!	setprngstates(0);
%!	sizeX = 5;
%!	sizeF = 4;
%!	vecX_secret = randn(sizeX,1);
%!	vecF_secret = randn(sizeF,1);
%!	matJ_secret = randn(sizeF,sizeX);
%!	funchF = @(x)( vecF_secret + matJ_secret*(x-vecX_secret) );
%!	vecX0 = randn(sizeX,1);
%!	%
%!	[ matJ0, ary3Kappa0, datOut ] = calcFDJKappa( vecX0, funchF );
%!	assert( isrealarray(matJ0,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
%!	assert( reldiff(matJ0,matJ_secret) < sqrt(eps) );
%!	%
%!	prm = [];
%!	[ matJ0, ary3Kappa0, datOut ] = calcFDJKappa( vecX0, funchF, prm );
%!	assert( isrealarray(matJ0,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
%!	assert( reldiff(matJ0,matJ_secret) < sqrt(eps) );
%!	%
%!	prm = [];
%!	prm.vecF0 = funchF( vecX0 );
%!	prm.vecEpsFD = sqrt(eps)*ones(sizeX,1);
%!	[ matJ0, ary3Kappa0, datOut ] = calcFDJKappa( vecX0, funchF, prm );
%!	assert( isrealarray(matJ0,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
%!	assert( reldiff(matJ0,matJ_secret) < sqrt(eps) );


%!function vecF = funcFQuad( vecX, vecX0, vecF0, matJ0, ary3Kappa0 )
%!	sizeX = size(vecX0,1);
%!	assert( isrealarray(vecX0,[sizeX,1]) );
%!	assert( isrealarray(vecX,[sizeX,1]) );
%!	sizeF = size(vecF0,1);
%!	assert( isrealarray(vecF0,[sizeF,1]) );
%!	assert( isrealarray(matJ0,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
%!	%
%!	vecD = vecX - vecX0;
%!	vecF = vecF0 + matJ0*vecD;
%!	for n=1:sizeF
%!		vecF(n) += 0.5*( vecD' * reshape( ary3Kappa0(n,:,:), [ sizeX, sizeX ] ) * vecD );
%!	endfor
%!endfunction

%!test
%!	setprngstates(0);
%!	sizeX = 2;
%!	sizeF = 2;
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
%!	funchF = @(x) funcFQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret );
%!	vecX0 = randn(sizeX,1);
%!	%
%!	[ matJ0, ary3Kappa0, datOut ] = calcFDJKappa( vecX0, funchF );
%!	assert( isrealarray(matJ0,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
%!	assert( reldiff(ary3Kappa0,ary3Kappa_secret) < eps^0.33 );
%!	%
%!	prm = [];
%!	[ matJ0, ary3Kappa0, datOut ] = calcFDJKappa( vecX0, funchF, prm );
%!	assert( isrealarray(matJ0,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
%!	assert( reldiff(ary3Kappa0,ary3Kappa_secret) < eps^0.33 );
%!	%
%!	prm = [];
%!	prm.vecF0 = funchF( vecX0 );
%!	prm.vecEpsFD = sqrt(sqrt(eps))*ones(sizeX,1);
%!	[ matJ0, ary3Kappa0, datOut ] = calcFDJKappa( vecX0, funchF, prm );
%!	assert( isrealarray(matJ0,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
%!	assert( reldiff(ary3Kappa0,ary3Kappa_secret) < eps^0.33 );


%!test
%!	setprngstates(0);
%!	sizeX = 2;
%!	sizeF = 2;
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
%!	funchF = @(x) funcFQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret );
%!	funchOmega = @(x) sumsq(funchF(x))/2.0;
%!	vecX0 = randn(sizeX,1);
%!	%
%!	[ matJ0, ary3Kappa0, datOut ] = calcFDJKappa( vecX0, funchF );
%!	assert( isrealarray(matJ0,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
%!	assert( reldiff(ary3Kappa0,ary3Kappa_secret) < eps^0.33 );
%!	%
%!	[ omega0, vecG0, matH0 ] = calcFDGH( vecX0, funchOmega );
%!	assert( reldiff(vecG0,datOut.vecG0) < eps^0.33 );
%!	assert( reldiff(matH0,datOut.matH0) < eps^0.33 );
%!	%echo__matH0 = matH0
%!	%echo__datOut_matH0 = datOut.matH0
