function [ funchXTrialOfP, curveGen_datOut ] = findZero_baseline__curveGen( vecX, vecF, matJ, ary3Kappa, funchF, curveGen_prm )
	curveGen_datOut = [];
	curveGen_datOut.fevalCount = 0;
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	matH = matJ'*matJ;
	[ matR, cholFlag ] = chol( matH );
	if ( 0 ~= cholFlag )
		error( "Hessian is not positive definite." );
	endif
	vecG = matJ'*vecF;
	matIX = eye(sizeX,sizeX);
	funchXTrialOfP = @(p)( vecX - ( p*matH + (1.0-p)*matIX ) \ (p*vecG) );
return;
endfunction
