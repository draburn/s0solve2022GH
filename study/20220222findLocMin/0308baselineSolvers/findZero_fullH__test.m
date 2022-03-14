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
	sizeX = 2;
	sizeF = sizeX;
	vecX_secret = randn(sizeX,1);
	vecF_secret = randn(sizeF,1);
	matJ_secret = randn(sizeF,sizeX);
	%
	addNullSpace = true;
	if ( addNullSpace )
		vecPhi_secret = randn(sizeX,1);
		vecPhi_secret /= norm(vecPhi_secret);
		matJ_secret = matJ_secret*( eye(sizeX,sizeX) - vecPhi_secret*(vecPhi_secret') );
	endif
	%
	sizeA = sizeX;
	for nf=1:sizeF
		matA = randn(sizeA,sizeX);
		matKappa = matA'*matA;
		makeKPositive = true;
		if (~makeKPositive)
			ary3Kappa_secret(nf,:,:) = matKappa;
		else
			ary3Kappa_secret(nf,:,:) = vecF_secret(nf)*matKappa;
		endif
	endfor
	%
	funchF = @(x) funcFQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret );
	vecX0 = randn(sizeX,1);
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
	powell_prm = [];
	powell_prm.stepCurveType = "powell";
	msg( __FILE__, __LINE__, "Calling findZero_fullH() POWELL..." );
	[ powell_vecXF, powell_datOut ] = findZero_fullH( vecX0, funchF, powell_prm );
	msg( __FILE__, __LINE__, "Finished findZero_fullH() POWELl." );
	%
	%
	newt_prm = [];
	newt_prm.stepCurveType = "newton";
	msg( __FILE__, __LINE__, "Calling findZero_fullH() NEWTON..." );
	[ newt_vecXF, newt_datOut ] = findZero_fullH( vecX0, funchF, newt_prm );
	msg( __FILE__, __LINE__, "Finished findZero_fullH() NEWTON." );
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
	omegaViz = min([ min(fsolveGnostic_datOut.omegaVals), min(lev_datOut.omegaVals) ]) -eps;
	%
	numFigs++; figure( numFigs );
	semilogy( ...
	  fsolveGnostic_datOut.fevalCountVals, fsolveGnostic_datOut.omegaVals-omegaViz, 'o-', 'markersize', 20, 'linewidth', 2, ...
	  lev_datOut.fevalCountVals, lev_datOut.omegaVals-omegaViz, 'x-', 'markersize', 20, 'linewidth', 2, ...
	  gradesc_datOut.fevalCountVals, gradesc_datOut.omegaVals-omegaViz, '+-', 'markersize', 20, 'linewidth', 2, ...
	  powell_datOut.fevalCountVals, powell_datOut.omegaVals-omegaViz, 's-', 'markersize', 20, 'linewidth', 2, ...
	  newt_datOut.fevalCountVals, newt_datOut.omegaVals-omegaViz, '^-', 'markersize', 20, 'linewidth', 2, ...
	  cauchy_datOut.fevalCountVals, cauchy_datOut.omegaVals-omegaViz, 'v-', 'markersize', 20, 'linewidth', 2 );
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
	  fsolveGnostic_datOut.fevalCountVals, fsolveGnostic_datOut.omegaVals-omegaViz, 'o-', 'markersize', 20, 'linewidth', 2, ...
	  lev_datOut.omegaVals-omegaViz, 'x-', 'markersize', 20, 'linewidth', 2, ...
	  gradesc_datOut.omegaVals-omegaViz, '+-', 'markersize', 20, 'linewidth', 2, ...
	  powell_datOut.omegaVals-omegaViz, 's-', 'markersize', 20, 'linewidth', 2, ...
	  newt_datOut.omegaVals-omegaViz, '^-', 'markersize', 20, 'linewidth', 2, ...
	  cauchy_datOut.omegaVals-omegaViz, 'v-', 'markersize', 20, 'linewidth', 2 );
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
