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
		assert( isposintscalar(nExactFit) );
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
	if ( isempty(s0) )
		s0 = xVals(nExactFit);
	end
	if ( isempty(p0) )
		p0 = 2.0;
	end
	if ( doChecks )
		assert( isrealscalar(s0) );
		assert( isrealscalar(p0) );
		assert( 0.0 < p0 );
	end
	%
	prm_calcAboutPt = mygetfield( prm, "prm_calcAboutPt", [] );
	[ bigF0, bigF1, rhoVals, omega, vecG, matH, matH2 ] = extFit__calcAboutPt( ...
	  s0, p0, xVals, fVals, nExactFit, wVals, prm_calcAboutPt );
	vecDeltaNewton = -matH\vecG;
	s1 = s0 + vecDeltaNewton(1);
	p1 = p0 + vecDeltaNewton(2);
	%
	%
	%
	genLevCurve = mygetfield( prm, "genLevCurve", true );
	if ( genLevCurve )
		tempCurve.vecG = vecG;
		tempCurve.matH = matH;
		tempCurve.matD = eye(2,2)*norm(diag(diag(matH)));
		tempCurve.numPts = 100;
		%
		tempCurve.matA = tempCurve.matH - tempCurve.matD;
		tempCurve.funchDeltaVec = @(lambda)( ...
		  -lambda * (( tempCurve.matD + lambda*tempCurve.matA ) \ vecG) );
		tempCurve.funchY = @(lambda) sqrt(sum( (tempCurve.funchDeltaVec(lambda)).^2 ));
		tempCurve.lambdaVals = daclinspace( 0.0, 1.0, tempCurve.numPts, tempCurve.funchY );
		for n=1:tempCurve.numPts
			tempCurve.vecDeltaVals(:,n) = tempCurve.funchDeltaVec(tempCurve.lambdaVals(n));
		end
		tempCurve.sVals = s0 + tempCurve.vecDeltaVals(1,:);
		tempCurve.pVals = p0 + tempCurve.vecDeltaVals(2,:);
		%
		levCurve = tempCurve;
		clear tempCurve;
	end
	%
	genLevMarqCurve = mygetfield( prm, "genLevMarqCurve", true );
	if ( genLevMarqCurve )
		tempCurve.vecG = vecG;
		tempCurve.matH = matH;
		tempCurve.matD = diag(diag(matH));
		tempCurve.numPts = 100;
		%
		tempCurve.matA = tempCurve.matH - tempCurve.matD;
		tempCurve.funchDeltaVec = @(lambda)( ...
		  -lambda * (( tempCurve.matD + lambda*tempCurve.matA ) \ vecG) );
		tempCurve.funchY = @(lambda) sqrt(sum( (tempCurve.funchDeltaVec(lambda)).^2 ));
		tempCurve.lambdaVals = daclinspace( 0.0, 1.0, tempCurve.numPts, tempCurve.funchY );
		for n=1:tempCurve.numPts
			tempCurve.vecDeltaVals(:,n) = tempCurve.funchDeltaVec(tempCurve.lambdaVals(n));
		end
		tempCurve.sVals = s0 + tempCurve.vecDeltaVals(1,:);
		tempCurve.pVals = p0 + tempCurve.vecDeltaVals(2,:);
		%
		levMarqCurve = tempCurve;
		clear tempCurve;
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
	sLo = min([ sLo, s0 ]);
	sHi = max([ sHi, s0 ]);
	pLo = min([ 1.0, p0 ]);
	pHi = max([ 10.0, p0 ]);
	sLo = mygetfield( prm, "sLo", sLo );
	sHi = mygetfield( prm, "sHi", sHi );
	pLo = mygetfield( prm, "pLo", pLo );
	pHi = mygetfield( prm, "pHi", pHi );
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
	zMesh = cap( zMesh, 0.0, (numContours+1.0)/(numContours-1.0)*max([ median(reshape(zMesh,1,[])), sqrt(omega) ]) );
	numFigs++; figure(numFigs);
	contourf( sVals_forMesh, pVals_forMesh, zMesh.^0.5, numContours );
	axis("square");
	%axis("equal");
	set( get(gcf,"children"), "ydir", "normal" );
	cmap = mycmap(numColors);
	colormap(cmap);
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( "sqrt(omega) vs s, p" );
	hold on;
	plot( ...
	  s0, p0, "s", "linewidth", 3, "markersize", 20, "color", [0.9,0.0,0.0], ...
	  s0, p0, "+", "linewidth", 3, "markersize", 20, "color", [0.9,0.0,0.0],
	  s1, p1, "s", "linewidth", 3, "markersize", 20, "color", [0.0,0.0,1.0], ...
	  s1, p1, "x", "linewidth", 3, "markersize", 20, "color", [0.0,0.0,1.0] );
	if ( genLevCurve )
		plot( levCurve.sVals, levCurve.pVals, "wo-" );
	end
	if ( genLevMarqCurve )
		plot( levMarqCurve.sVals, levMarqCurve.pVals, "ko-" );
	end
	hold off;
return;
end
