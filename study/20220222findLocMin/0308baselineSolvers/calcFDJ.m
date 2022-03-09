% Function...
%  Calculates the Jacobian of a function via finite-differencing.
 
function [ matJ, datOut ] = calcFDJ( vecX0, funchF, prm=[] )
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
	fdOrder = mygetfield( prm, "fdOrder", 1 );
	assert( isrealscalar(fdOrder) );
	%
	vecEpsFD = mygetfield( prm, "vecEpsFD", [] );
	if ( isempty(vecEpsFD) )
		switch (fdOrder)
		case 1
			epsFDRel = eps^0.5;
		case 2
			epsFDRel = eps^0.33;
		otherwise
			error( "Invalid value of fdOrder." );
		endswitch
		epsFDRel = mygetfield( prm, "epsFDRel", epsFDRel );
		assert( isrealscalar(epsFDRel) );
		assert( 0.0 < epsFDRel );
		vecEpsFD = epsFDRel * ( 1.0 + abs(vecX0) );
	endif
	assert( isrealarray(vecEpsFD,[sizeX,1]) );
	assert( sum(double(0.0<=vecEpsFD)) == sizeX );
	%
	%
	%
	matJ = zeros(sizeF,sizeX);
	switch (fdOrder)
	case 1
		for n = 1 : sizeX
			vecXP = vecX0;
			vecXP(n) += vecEpsFD(n);
			vecFP = funchF(vecXP); fevalCount++;
			matJ(:,n) = (vecFP-vecF0)/vecEpsFD(n);
		endfor
	case 2
		for n = 1 : sizeX
			vecXP = vecX0;
			vecXM = vecX0;
			vecXP(n) += vecEpsFD(n);
			vecXM(n) -= vecEpsFD(n);
			vecFP = funchF(vecXP); fevalCount++;
			vecFM = funchF(vecXM); fevalCount++;
			matJ(:,n) = (vecFP-vecFM)/(2.0*vecEpsFD(n));
		endfor
	otherwise
		error( "Invalid value of fdOrder." );
	endswitch
	%
	if ( nargout >= 2 )
		datOut.fevalCount = fevalCount;
		datOut.fdOrder = fdOrder;
		datOut.vecEpsFD = vecEpsFD;
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
%!	[ matJ, datOut ] = calcFDJ( vecX0, funchF );
%!	fevalCount = datOut.fevalCount
%!	assert( reldiff(matJ,matJ_secret) < sqrt(eps) );
%!	assert( fevalCount == 1+sizeX );
%!	%
%!	prm = [];
%!	[ matJ, datOut ] = calcFDJ( vecX0, funchF, prm );
%!	fevalCount = datOut.fevalCount
%!	assert( reldiff(matJ,matJ_secret) < sqrt(eps) );
%!	assert( fevalCount == 1+sizeX );
%!	%
%!	prm = [];
%!	prm.vecF0 = funchF( vecX0 );
%!	prm.fdOrder = 2;
%!	prm.vecEpsFD = sqrt(eps)*ones(sizeX,1);
%!	[ matJ, datOut ] = calcFDJ( vecX0, funchF, prm );
%!	fevalCount = datOut.fevalCount
%!	assert( reldiff(matJ,matJ_secret) < sqrt(eps) );
%!	assert( fevalCount == 2*sizeX );
