% Function...

function [ matSC, vecGSC, matHSC, datOut ] = findZero_baseline__scale( vecG, matH, prm=[] )
	sizeX = size(vecG,1);
	matSC = mygetfield( prm, "matSC", [] );
	if ( isempty(matSC) )
		stepScalingType = mygetfield( prm, "stepScalingType", "" );
		switch ( tolower(stepScalingType) )
		case { "", "none", "identity" }
			matSC = eye(sizeX,sizeX);
		case { "marquardt" }
			hNorm = sqrt( sum(sumsq(matH))/sizeX );
			matSC = diag( 1.0./sqrt( abs(diag(matH)) + eps*hNorm ) );
		otherwise
			error( "Invalid value of stepScalingType." );
		endswitch
	endif
	assert( isrealarray(matSC,[sizeX,sizeX]) );
	vecGSC = matSC'*vecG;
	matHSC = matSC'*matH*matSC;
return;
endfunction


%!test
%!	setprngstates(0);
%!	numFigs = 0;
%!	%
%!	sizeX = 2;
%!	sizeF = 3;
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX);
%!	omega = sumsq(vecF)/2.0;
%!	vecG = matJ'*vecF
%!	matH = matJ'*matJ
%!	%
%!	[ matSC, vecGSC, matHSC ] = findZero_baseline__scale( vecG, matH );
%!	%
%!	prm = []; prm.matSC = randn(sizeX,sizeX)
%!	[ matSC, vecGSC, matHSC ] = findZero_baseline__scale( vecG, matH, prm )
%!	%
%!	prm = []; prm.stepScalingType = "none"
%!	[ matSC, vecGSC, matHSC ] = findZero_baseline__scale( vecG, matH, prm )
%!	%
%!	prm = []; prm.stepScalingType = "marquardt"
%!	[ matSC, vecGSC, matHSC ] = findZero_baseline__scale( vecG, matH, prm )
