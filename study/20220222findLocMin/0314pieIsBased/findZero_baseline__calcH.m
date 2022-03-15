% Function...
 
function [ matH, datOut ] = findZero_baseline__calcH( vecX, vecF, matJ, funchF, prm=[] )
	matJTJ = matJ'*matJ;
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	fevalCount = 0;
	cvtkv_prm = mygetfield( prm, "cvtkv_prm", [] );
	kMode = mygetfield( prm, "kMode", "none" );
	switch ( tolower(kMode) )
	case { "", "none", "jtj" }
		matV = [];
	case { "patch", "phi" }
		cvtkv_prm.includeOffDiagonalTerms = false;
		patchTol = mygetfield( prm, "patchTol", eps^0.5 );
		[ matPsi_jtj, matLambda_jtj ] = eig( matJTJ );
		vecLambda_jtj = diag(matLambda_jtj);
		m = ( vecLambda_jtj < patchTol*max(abs(vecLambda_jtj)) );
		matV = matPsi_jtj(:,m);
	case { "dsb", "diagonal standard basis" }
		cvtkv_prm.includeOffDiagonalTerms = false;
		matV = eye(sizeX,sizeX);
	case { "deb", "diagonal eigenbasis" }
		cvtkv_prm.includeOffDiagonalTerms = false;
		[ matPsi_jtj, matLambda_jtj ] = eig( matJTJ );
		matV = matPsi_jtj;
	case { "full" }
		cvtkv_prm.includeOffDiagonalTerms = true;
		matV = eye(sizeX,sizeX);
	otherwise
		error( "Invalid value of kMode." );
	endswitch
	sizeV = size(matV,2);
	if ( 0 == sizeV )
		matH = matJTJ;
		ary3KappaA = zeros(sizeF,sizeX,sizeX);
	else
		cvtkv_prm.calcSBKappaA = true;
		[ matVTKV, cvtkv_datOut ] = calcFDVTKV( vecX, vecF, funchF, matV, cvtkv_prm );
		fevalCount += cvtkv_datOut.fevalCount;
		matH = matJTJ + (matV * matVTKV * (matV'));
		ary3KappaA = cvtkv_datOut.ary3SBKappaA;
	endif
	if ( nargout >= 2 )
		datOut.fevalCount = fevalCount;
		datOut.ary3KappaA = ary3KappaA;
	endif
return;
endfunction


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
%!	funchF = @(x) funcFQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret);
%!	vecX = randn(sizeX,1);
%!	vecF = funchF(vecX);
%!	%
%!	cfdj_prm = [];
%!	cfdj_prm.fdOrder = 2;
%!	matJ = calcFDJ( vecX, funchF, cfdj_prm );
%!	%
%!	matH = findZero_baseline__calcH( vecX, vecF, matJ, funchF )
%!	rcond(matH)
%!	%
%!	ch_prm = []; ch_prm.kMode = "patch";
%!	matH = findZero_baseline__calcH( vecX, vecF, matJ, funchF, ch_prm )
%!	rcond(matH)
%!	%
%!	ch_prm = []; ch_prm.kMode = "diagonal standard basis";
%!	matH = findZero_baseline__calcH( vecX, vecF, matJ, funchF, ch_prm )
%!	rcond(matH)
%!	%
%!	ch_prm = []; ch_prm.kMode = "diagonal eigenbasis";
%!	matH = findZero_baseline__calcH( vecX, vecF, matJ, funchF, ch_prm )
%!	rcond(matH)
%!	%
%!	ch_prm = []; ch_prm.kMode = "full";
%!	matH = findZero_baseline__calcH( vecX, vecF, matJ, funchF, ch_prm )
%!	rcond(matH)
