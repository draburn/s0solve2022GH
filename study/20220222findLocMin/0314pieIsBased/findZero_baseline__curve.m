% Function...

function [ funchDeltaOfP, datOut ] = findZero_baseline__curve( omega, vecGSC, matHSC, matHSCRegu, matSC, prm=[] )
	vecDeltaNewton = matSC * ( matHSCRegu \ (-vecGSC) );
	pCauchy = calcLinishRootOfQuad( 0.5*(vecGSC'*matHSC*vecGSC), -sumsq(vecGSC), omega );
	assert( pCauchy > 0.0 );
	vecDeltaCauchy = pCauchy*matSC*(-vecGSC);
	sizeX = size(vecGSC,1);
	stepCurveType = mygetfield( prm, "stepCurveType", "" );
	switch ( tolower(stepCurveType) )
	case { "newton" }
		funchDeltaOfP = @(p) ( p * vecDeltaNewton );
	case { "", "levenberg" }
		funchDeltaOfP = @(p) matSC * (( p * matHSCRegu + (1.0-p)*eye(sizeX,sizeX) ) \ (-p*vecGSC));
	case { "powell", "dog leg" }
		funchDeltaOfP = @(p) ( 2.0*p*vecDeltaCauchy + ...
		  (p>0.5) * ( (2.0*p-1.0)*vecDeltaNewton + 4.0*(0.5-p)*vecDeltaCauchy ) );
	case { "gradesc", "gradient descent curve" }
		% I suspect we could get a faster run time by doing an ODE solve then interpolating, but, POITROME.
		[ matPsi_hscRegu, matLambda_hscRegu ] = eig( matHSCRegu );
		vecLambda_hscRegu = diag(matLambda_hscRegu);
		vecLIPNG = matLambda_hscRegu \ ( matPsi_hscRegu' * (-vecGSC) );
		matSP = matSC * matPsi_hscRegu;
		funchDeltaOfP = @(p) ( matSP * (diag( 1.0 - (1.0-p).^vecLambda_hscRegu ) * vecLIPNG) );
	case { "cauchy", "gradient descent segment" }
		funchDeltaOfP = @(p) ( p * vecDeltaCauchy );
	otherwise
		error( "Invalid value of stepCurveType." );
	endswitch
	if ( nargout >= 2 )
		datOut.vecDeltaNewton = vecDeltaNewton;
		datOut.pCauchy = pCauchy;
		datOut.vecDeltaCauchy = vecDeltaCauchy;
	endif
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
%!	vecGSC = matJ'*vecF;
%!	matHSC = matJ'*matJ;
%!	matHSCRegu = findZero_baseline__regularize( matHSC );
%!	matSC = [ 1.0, 0.0; 0.0, 1.0 ];
%!	%
%!	prm = []; prm.stepCurveType = "levenberg";
%!	funchDeltaOfP_lev = findZero_baseline__curve( omega, vecGSC, matHSC, matHSCRegu, matSC, prm );
%!	%
%!	prm = []; prm.stepCurveType = "gradesc";
%!	funchDeltaOfP_gradesc = findZero_baseline__curve( omega, vecGSC, matHSC, matHSCRegu, matSC, prm );
%!	%
%!	prm = []; prm.stepCurveType = "powell";
%!	funchDeltaOfP_powell = findZero_baseline__curve( omega, vecGSC, matHSC, matHSCRegu, matSC, prm );
%!	%
%!	prm = []; prm.stepCurveType = "newton";
%!	funchDeltaOfP_newton = findZero_baseline__curve( omega, vecGSC, matHSC, matHSCRegu, matSC, prm );
%!	%
%!	prm = []; prm.stepCurveType = "cauchy";
%!	funchDeltaOfP_cauchy = findZero_baseline__curve( omega, vecGSC, matHSC, matHSCRegu, matSC, prm );
%!	%
%!	numPts = 101;
%!	pPts = linspace( 0.0, 1.0, numPts );
%!	for n=1:numPts
%!		vecDeltaPts_lev(:,n) = funchDeltaOfP_lev(pPts(n));
%!		vecDeltaPts_gradesc(:,n) = funchDeltaOfP_gradesc(pPts(n));
%!		vecDeltaPts_powell(:,n) = funchDeltaOfP_powell(pPts(n));
%!		vecDeltaPts_newton(:,n) = funchDeltaOfP_newton(pPts(n));
%!		vecDeltaPts_cauchy(:,n) = funchDeltaOfP_cauchy(pPts(n));
%!	endfor
%!	%
%!	numFigs = 0;
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  0, 0, 'p', 'linewidth', 3, 'markersize', 25, ...
%!	  vecDeltaPts_lev(1,:), vecDeltaPts_lev(2,:), 'o-', ...
%!	  vecDeltaPts_gradesc(1,:), vecDeltaPts_gradesc(2,:), 'o-', ...
%!	  vecDeltaPts_powell(1,:), vecDeltaPts_powell(2,:), 'o-', ...
%!	  vecDeltaPts_newton(1,:), vecDeltaPts_newton(2,:), 'o-', ...
%!	  vecDeltaPts_cauchy(1,:), vecDeltaPts_cauchy(2,:), 'o-', ...
%!	  vecDeltaPts_newton(1,end), vecDeltaPts_newton(2,end), 'x', 'linewidth', 3, 'markersize', 25, ...
%!	  vecDeltaPts_cauchy(1,end), vecDeltaPts_cauchy(2,end), '+', 'linewidth', 3, 'markersize', 25 );
%!	grid on;
