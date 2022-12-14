	clear
	numFigs = 0;
	setprngstates(0);
	%
	function vecF = funcFQuad( vecX, vecX0, vecF0, matJ0, ary3Kappa0 )
		sizeX = size(vecX0,1);
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecX,[sizeX,1]) );
		sizeF = size(vecF0,1);
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
		%
		vecD = vecX - vecX0;
		vecF = vecF0 + matJ0*vecD;
		for n=1:sizeF
			vecF(n) += 0.5*( vecD' * reshape( ary3Kappa0(n,:,:), [ sizeX, sizeX ] ) * vecD );
		endfor
	endfunction
	%
	caseNum = 200;
	msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
	switch (caseNum)
	case 100
		sizeX = 3;
		sizeF = sizeX;
		vecX_secret = randn(sizeX,1);
		vecF_secret = zeros(sizeF,1);
		matJ_secret = randn(sizeF,sizeX);
		for nf=1:sizeF
			matKappa = randn(sizeX,sizeX);
			matKappa = ( matKappa' + matKappa ) / 2.0;
			ary3Kappa_secret(nf,:,:) = matKappa;
		endfor
		funchF = @(x) funcFQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret );
		vecX0 = randn(sizeX,1);
	case 200
		sizeX = 5;
		sizeF = sizeX;
		vecX_secret = randn(sizeX,1);
		matJ_secret = randn(sizeF,sizeX);
		matA0 = randn(sizeF,sizeX);
		matA1 = randn(sizeX,sizeX);
		matA2 = randn(sizeX,sizeX);
		matB0 = randn(sizeF,sizeX);
		matB1 = randn(sizeX,sizeX);
		matB2 = randn(sizeX,sizeX);
		matB3 = randn(sizeX,sizeX);
		y = @(x)( x - vecX_secret );
		funchF = @(x)( matJ_secret*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
		vecX0 = randn(sizeX,1);
	otherwise
		error( "Invalid case." );
	endswitch
	%
	%
	%
	%
	%
	lev_prm = [];
	msg( __FILE__, __LINE__, "Calling findZero_fullH() LEVENBERG..." );
	[ lev_vecXF, lev_datOut ] = findZero_fullH( vecX0, funchF, lev_prm );
	msg( __FILE__, __LINE__, "Finished findZero_fullH() LEVENBERG." );
	%
	%
	gradesc_prm = [];
	gradesc_prm.stepCurveType = "gradient descent curve";
	msg( __FILE__, __LINE__, "Calling findZero_fullH() GRADIENT DESCENT CURVE..." );
	[ gradesc_vecXF, gradesc_datOut ] = findZero_fullH( vecX0, funchF, gradesc_prm );
	msg( __FILE__, __LINE__, "Finished findZero_fullH() GRADIENT DESCENT CURVE." );
	%
	%
	newt_prm = [];
	newt_prm.stepCurveType = "newton";
	msg( __FILE__, __LINE__, "Calling findZero_fullH() NEWTON..." );
	[ newt_vecXF, newt_datOut ] = findZero_fullH( vecX0, funchF, newt_prm );
	msg( __FILE__, __LINE__, "Finished findZero_fullH() NEWTON." );
	%
	%
	powell_prm = [];
	powell_prm.stepCurveType = "powell";
	msg( __FILE__, __LINE__, "Calling findZero_fullH() POWELL..." );
	[ powell_vecXF, powell_datOut ] = findZero_fullH( vecX0, funchF, powell_prm );
	msg( __FILE__, __LINE__, "Finished findZero_fullH() POWELl." );
	%
	%
	cauchy_prm = [];
	cauchy_prm.stepCurveType = "cauchy";
	msg( __FILE__, __LINE__, "Calling findZero_fullH() CAUCHY..." );
	[ cauchy_vecXF, cauchy_datOut ] = findZero_fullH( vecX0, funchF, cauchy_prm );
	msg( __FILE__, __LINE__, "Finished findZero_fullH() CAUCHY." );
	%
	%
	msg( __FILE__, __LINE__, "Calling findZero_fsolveGnostic()..." );
	fsolveGnostic_prm = [];
	[ fsolveGnostic_vecXF, fsolveGnostic_datOut ] = findZero_fsolveGnostic( vecX0, funchF, fsolveGnostic_prm );
	msg( __FILE__, __LINE__, "Finished findZero_fsolveGnostic()." );
	%
	%
	omegaViz = min([ ...
	  min(lev_datOut.omegaVals), ...
	  min(gradesc_datOut.omegaVals), ...
	  min(powell_datOut.omegaVals), ...
	  min(newt_datOut.omegaVals), ...
	  min(cauchy_datOut.omegaVals), ...
	  min(fsolveGnostic_datOut.omegaVals) ]);
	%
	numFigs++; figure( numFigs );
	semilogy( ...
	  fsolveGnostic_datOut.fevalCountVals, fsolveGnostic_datOut.omegaVals-omegaViz+eps^2, 'o-', 'markersize', 20, 'linewidth', 2, ...
	  lev_datOut.fevalCountVals, lev_datOut.omegaVals-omegaViz+eps^2, 'x-', 'markersize', 20, 'linewidth', 2, ...
	  gradesc_datOut.fevalCountVals, gradesc_datOut.omegaVals-omegaViz+eps^2, '+-', 'markersize', 20, 'linewidth', 2, ...
	  powell_datOut.fevalCountVals, powell_datOut.omegaVals-omegaViz+eps^2, 's-', 'markersize', 20, 'linewidth', 2, ...
	  newt_datOut.fevalCountVals, newt_datOut.omegaVals-omegaViz+eps^2, '^-', 'markersize', 20, 'linewidth', 2, ...
	  cauchy_datOut.fevalCountVals, cauchy_datOut.omegaVals-omegaViz+eps^2, 'v-', 'markersize', 20, 'linewidth', 2 );
	grid on;
	ylabel( "omega-omegaViz" );
	xlabel( "feval count" );
	legend( ...
	  "fsolve gnostic", ...
	  "levenberg", ...
	  "gradesc", ...
	  "powell", ...
	  "newton", ...
	  "cauchy", ...
	  "location", "northEast" );
	%
	numFigs++; figure( numFigs );
	semilogy( ...
	  fsolveGnostic_datOut.fevalCountVals, fsolveGnostic_datOut.omegaVals-omegaViz+eps^2, 'o-', 'markersize', 20, 'linewidth', 2, ...
	  lev_datOut.omegaVals-omegaViz+eps^2, 'x-', 'markersize', 20, 'linewidth', 2, ...
	  gradesc_datOut.omegaVals-omegaViz+eps^2, '+-', 'markersize', 20, 'linewidth', 2, ...
	  powell_datOut.omegaVals-omegaViz+eps^2, 's-', 'markersize', 20, 'linewidth', 2, ...
	  newt_datOut.omegaVals-omegaViz+eps^2, '^-', 'markersize', 20, 'linewidth', 2, ...
	  cauchy_datOut.omegaVals-omegaViz+eps^2, 'v-', 'markersize', 20, 'linewidth', 2 );
	grid on;
	ylabel( "omega-omegaViz" );
	xlabel( "iter count or feval count" );
	legend( ...
	  "fsolve gnostic", ...
	  "levenberg", ...
	  "gradesc", ...
	  "powell", ...
	  "newton", ...
	  "cauchy", ...
	  "location", "northEast" );
	%
	msg( __FILE__, __LINE__, "Please see figure(s)." );
