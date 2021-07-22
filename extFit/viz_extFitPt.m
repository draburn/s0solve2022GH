function viz_extFitPt( xVals, fVals, nExactFit, s0=[], p0=[], wVals=[], prm=[] )
	commondefs;
	thisFile = "viz_extFitPt";
	doChecks = mygetfield( prm, "doChecks", true );
	numFigs = mygetfield( prm, "numFigs0", 0 );
	%
	if ( isempty(wVals) )
		wVals = ones(size(xVals));
	end
	if ( doChecks )
		numPts = size(xVals,2);
		assert( numPts >= 3 );
		assert( isrealarray(xVals,[1,numPts]) );
		xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
		assert( xValsAreStrictlyIncreasing );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isrealscalar(nExactFit) );
		assert( fleq(nExactFit,round(nExactFit)) );
		assert( 1 <= nExactFit );
		assert( nExactFit <= numPts );
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	%
	%
	if ( ~isempty(s0) || ~isempty(p0) )
		if ( doChecks )
			assert( isrealscalar(s0) );
			assert( isrealscalar(p0) );
			assert( 0.0 < p0 );
		end
		prm_calcAboutPt = mygetfield( prm, "prm_calcAboutPt", [] );
		[ bigF0, bigF1, rhoVals, omega, vecG, matH, matH2 ] = extFit__calcAboutPt( ...
		  s0, p0, xVals, fVals, nExactFit, wVals, prm_calcAboutPt );
		matD = diag(abs(diag(matH)));
		%
		vecDelta = -matH\vecG;
		s1 = s0 + vecDelta(1);
		p1 = p0 + vecDelta(2);
		%
		numCurvePts = 100;
		matA = matH - matD;
		funchY_forCurve = @(x)(norm(x*( (matD+x*matA)\vecG )));
		lambdaVals_forCurve = daclinspace( 0.0, 1.0, numCurvePts, funchY_forCurve );
		for n = 1 : numCurvePts
			lambda = lambdaVals_forCurve(n);
			if ( eps >= lambda )
				vecDelta = -lambda*(matD\vecG);
			else
				mu = (1.0/lambda)-1.0;
				vecDelta = -( matH + mu*matD ) \ vecG;
			end
			sVals_forCurve(n) = s0 + vecDelta(1);
			pVals_forCurve(n) = p0 + vecDelta(2);
		end
	end
	%
	%
	%
	if ( 1 == nExactFit )
		sLo = xVals(1) - 1.0*(xVals(numPts)-xVals(1));
		sHi = xVals(2);
	elseif ( numPts == nExactFit )
		sLo = xVals(numPts-1);
		sHi = xVals(numPts) + 1.0*(xVals(numPts)-xVals(1));
	else
		sLo = xVals(nExactFit-1);
		sHi = xVals(nExactFit+1);
	end
	sLo = mygetfield( prm, "sLo", sLo );
	sHi = mygetfield( prm, "sHi", sHi );
	pLo = mygetfield( prm, "pLo", 1.0 );
	pHi = mygetfield( prm, "pHi", 10.0 );
	sSize = mygetfield( prm, "sSize", 200 );
	pSize = mygetfield( prm, "pSize", 201 );
	if ( doChecks )
		assert( isrealscalar(sLo) );
		assert( isrealscalar(sHi) );
		assert( isrealscalar(pLo) );
		assert( isrealscalar(pHi) );
		assert( isposintscalar(sSize) );
		assert( isposintscalar(pSize) );
		assert( sLo < sHi );
		assert( 0.0 < pLo );
		assert( pLo < pHi );
		assert( 3 <= sSize );
		assert( 3 <= pSize );
	end
	%
	sVals_forMesh = linspace( sLo, sHi, sSize );
	pVals_forMesh = linspace( pLo, pHi, pSize );
	[ sMesh, pMesh ] = meshgrid( sVals_forMesh, pVals_forMesh );
	size1 = size(sMesh,1);
	size2 = size(sMesh,2);
	%
	prm_calcMesh = mygetfield( prm, "prm_calcMesh", [] );
	[ bigF0Mesh, bigF1Mesh, omegaMesh ] = extFit__calcMesh( ...
	  sMesh, pMesh, xVals, fVals, nExactFit, wVals, prm_calcMesh );
	%
	%
	%
	numContours = mygetfield( prm, "numContours", 20 );
	numColors = mygetfield( prm, "numColors", 21 );
	if ( doChecks )
		assert( isposintscalar(numContours) );
		assert( isposintscalar(numColors) );
		assert( 3 <= numContours );
		assert( 3 <= numColors );
	end
	%
	%
	%
	zMesh = sqrt(omegaMesh);
	numFigs++; figure(numFigs);
	contourf( sVals_forMesh, pVals_forMesh, zMesh.^0.5, numContours );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	cmap = mycmap(numColors);
	colormap(cmap);
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( "sqrt(omega) vs s, p" );
	if ( ~isempty(s0) || ~isempty(p0) )
		hold on; % Continue omega contour graph.
		plot( ...
		  s0, p0, "s", "linewidth", 3, "markersize", 20, "color", [0.9,0.0,0.0], ...
		  s0, p0, "+", "linewidth", 3, "markersize", 20, "color", [0.9,0.0,0.0],
		  s1, p1, "s", "linewidth", 3, "markersize", 20, "color", [0.0,0.0,1.0], ...
		  s1, p1, "x", "linewidth", 3, "markersize", 20, "color", [0.0,0.0,1.0] );
		plot( sVals_forCurve, pVals_forCurve, "ko-" );
		hold off;
	end
return;
end
