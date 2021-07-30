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
	[ rhoVals, bigF0, bigF1, omega, vecG, matH ] = extFit__calcAboutPt( ...
	  s0, p0, xVals, fVals, nExactFit, wVals, prm_calcAboutPt );
	vecDeltaNewton = -matH\vecG;
	s1 = s0 + vecDeltaNewton(1);
	p1 = p0 + vecDeltaNewton(2);
	rvecHTHDiag = sum(matH.^2,1);
	hScale = sqrt(sum(rvecHTHDiag)/2.0);
	%
	%
	%
	prm_calcMesh = mygetfield( prm, "prm_calcMesh", [] );
	prm_daclinspace = [];
	prm_daclinspace.coeffMin = 0.9;
	prm_daclinspace.coeffMax = 1.1;
	prm_daclinspace.numIterLimit = 100;
	genLevCurve = mygetfield( prm, "genLevCurve", true );
	if ( genLevCurve )
		tempCurve.omega0 = omega;
		tempCurve.vecG = vecG;
		tempCurve.matH = matH;
		tempCurve.matD = eye(2,2)*norm(diag(diag(matH)));
		tempCurve.numPts = 1000;
		%
		tempCurve.matA = tempCurve.matH - tempCurve.matD;
		tempCurve.funchDeltaVec = @(lambda)( ...
		  -lambda * (( tempCurve.matD + lambda*tempCurve.matA ) \ vecG) );
		if (0)
			tempCurve.lambdaVals = daclinspace( 0.0, 1.0, tempCurve.numPts, ...
			  tempCurve.funchDeltaVec, prm_daclinspace );
		else
			tempCurve.muVals = [ -1.0, eps.^linspace(-0.4,0.4,tempCurve.numPts-2), 0.0 ];
			tempCurve.lambdaVals = [ 0.0, 1.0./(1.0+tempCurve.muVals(2:end)) ];
		end
		%
		for n=1:tempCurve.numPts
			tempCurve.vecDeltaVals(:,n) = tempCurve.funchDeltaVec(tempCurve.lambdaVals(n));
			tempCurve.omegaModelVals(:,n) = tempCurve.omega0 ...
			 + tempCurve.vecDeltaVals(:,n)'*tempCurve.vecG ...
			 + 0.5*tempCurve.vecDeltaVals(:,n)'*tempCurve.matH*+ tempCurve.vecDeltaVals(:,n);
		end
		tempCurve.sVals = s0 + tempCurve.vecDeltaVals(1,:);
		tempCurve.pVals = p0 + tempCurve.vecDeltaVals(2,:);
		tempCurve.dacVals(1) = 0.0;
		for n=2:tempCurve.numPts
			tempCurve.dacVals(n) = tempCurve.dacVals(n-1) + ...
			  norm( tempCurve.vecDeltaVals(:,n) - tempCurve.vecDeltaVals(:,n-1) );
		end
		tempCurve.deltaNormVals = sqrt(sum(tempCurve.vecDeltaVals.^2,1));
		[ tempCurve.bigF0Vals, tempCurve.bigF1Vals, tempCurve.omegaVals ] = extFit__calcMesh( ...
		  tempCurve.sVals, tempCurve.pVals, xVals, fVals, nExactFit, wVals, prm_calcMesh );
		%tempCurve.alphaVals = -(tempCurve.matD\tempCurve.vecG)' * tempCurve.vecDeltaVals ...
		%  ./ ( norm(tempCurve.matD\tempCurve.vecG) * tempCurve.deltaNormVals );
		tempCurve.vecDDeltaVals = tempCurve.vecDeltaVals(:,2:end) - tempCurve.vecDeltaVals(:,1:end-1);
		tempCurve.alphaVals = (-(tempCurve.matD\tempCurve.vecG)' * tempCurve.vecDDeltaVals) ...
		 ./ ( norm(tempCurve.matD\tempCurve.vecG) * sqrt(sum(tempCurve.vecDDeltaVals.^2,1)) );
		%
		levCurve = tempCurve;
		clear tempCurve;
	end
	%
	genLevMarqCurve = mygetfield( prm, "genLevMarqCurve", true );
	if ( genLevMarqCurve )
		tempCurve.omega0 = omega;
		tempCurve.vecG = vecG;
		tempCurve.matH = matH;
		tempCurve.matD = diag(diag(matH));
		tempCurve.numPts = 1000;
		%
		tempCurve.matA = tempCurve.matH - tempCurve.matD;
		tempCurve.funchDeltaVec = @(lambda)( ...
		  -lambda * (( tempCurve.matD + lambda*tempCurve.matA ) \ vecG) );
		if (0)
			tempCurve.lambdaVals = daclinspace( 0.0, 1.0, tempCurve.numPts, ...
			  tempCurve.funchDeltaVec, prm_daclinspace );
		else
			tempCurve.muVals = [ -1.0, eps.^linspace(-0.4,0.4,tempCurve.numPts-2), 0.0 ];
			tempCurve.lambdaVals = [ 0.0, 1.0./(1.0+tempCurve.muVals(2:end)) ];
		end
		%
		for n=1:tempCurve.numPts
			tempCurve.vecDeltaVals(:,n) = tempCurve.funchDeltaVec(tempCurve.lambdaVals(n));
			tempCurve.omegaModelVals(:,n) = tempCurve.omega0 ...
			 + tempCurve.vecDeltaVals(:,n)'*tempCurve.vecG ...
			 + 0.5*tempCurve.vecDeltaVals(:,n)'*tempCurve.matH*+ tempCurve.vecDeltaVals(:,n);
		end
		tempCurve.sVals = s0 + tempCurve.vecDeltaVals(1,:);
		tempCurve.pVals = p0 + tempCurve.vecDeltaVals(2,:);
		tempCurve.dacVals(1) = 0.0;
		for n=2:tempCurve.numPts
			tempCurve.dacVals(n) = tempCurve.dacVals(n-1) + ...
			  norm( tempCurve.vecDeltaVals(:,n) - tempCurve.vecDeltaVals(:,n-1) );
		end
		tempCurve.deltaNormVals = sqrt(sum(tempCurve.vecDeltaVals.^2,1));
		[ tempCurve.bigF0Vals, tempCurve.bigF1Vals, tempCurve.omegaVals ] = extFit__calcMesh( ...
		  tempCurve.sVals, tempCurve.pVals, xVals, fVals, nExactFit, wVals, prm_calcMesh );
		%tempCurve.alphaVals = -(tempCurve.matD\tempCurve.vecG)' * tempCurve.vecDeltaVals ...
		%  ./ ( norm(tempCurve.matD\tempCurve.vecG) * tempCurve.deltaNormVals );
		tempCurve.vecDDeltaVals = tempCurve.vecDeltaVals(:,2:end) - tempCurve.vecDeltaVals(:,1:end-1);
		tempCurve.alphaVals = (-(tempCurve.matD\tempCurve.vecG)' * tempCurve.vecDDeltaVals) ...
		 ./ ( norm(tempCurve.matD\tempCurve.vecG) * sqrt(sum(tempCurve.vecDDeltaVals.^2,1)) );
		%
		levMarqCurve = tempCurve;
		clear tempCurve;
	end
	%
	reguCoeff = mygetfield( prm, "reguCoeff", eps );
	if ( 1 )
		tempCurve.omega0 = omega;
		tempCurve.vecG = vecG;
		tempCurve.matH = matH;
		matDMod =  sqrt( diag(rvecHTHDiag) + (hScale^2)*reguCoeff*eye(2,2) );
		tempCurve.matD = matDMod;
		%%%tempCurve.matD = diag(abs(diag(matH))) + hScale*reguCoeff*eye(2,2);
		tempCurve.numPts = 1000;
		%
		tempCurve.matA = tempCurve.matH - tempCurve.matD;
		tempCurve.funchDeltaVec = @(lambda)( ...
		  -lambda * (( tempCurve.matD + lambda*tempCurve.matA ) \ vecG) );
		if (0)
			tempCurve.lambdaVals = daclinspace( 0.0, 1.0, tempCurve.numPts, ...
			  tempCurve.funchDeltaVec, prm_daclinspace );
		else
			tempCurve.muVals = [ -1.0, eps.^linspace(-0.4,0.4,tempCurve.numPts-2), 0.0 ];
			tempCurve.lambdaVals = [ 0.0, 1.0./(1.0+tempCurve.muVals(2:end)) ];
		end
		%
		for n=1:tempCurve.numPts
			tempCurve.vecDeltaVals(:,n) = tempCurve.funchDeltaVec(tempCurve.lambdaVals(n));
			tempCurve.omegaModelVals(:,n) = tempCurve.omega0 ...
			 + tempCurve.vecDeltaVals(:,n)'*tempCurve.vecG ...
			 + 0.5*tempCurve.vecDeltaVals(:,n)'*tempCurve.matH*+ tempCurve.vecDeltaVals(:,n);
		end
		tempCurve.sVals = s0 + tempCurve.vecDeltaVals(1,:);
		tempCurve.pVals = p0 + tempCurve.vecDeltaVals(2,:);
		tempCurve.dacVals(1) = 0.0;
		for n=2:tempCurve.numPts
			tempCurve.dacVals(n) = tempCurve.dacVals(n-1) + ...
			  norm( tempCurve.vecDeltaVals(:,n) - tempCurve.vecDeltaVals(:,n-1) );
		end
		tempCurve.deltaNormVals = sqrt(sum(tempCurve.vecDeltaVals.^2,1));
		[ tempCurve.bigF0Vals, tempCurve.bigF1Vals, tempCurve.omegaVals ] = extFit__calcMesh( ...
		  tempCurve.sVals, tempCurve.pVals, xVals, fVals, nExactFit, wVals, prm_calcMesh );
		%tempCurve.alphaVals = -(tempCurve.matD\tempCurve.vecG)' * tempCurve.vecDeltaVals ...
		%  ./ ( norm(tempCurve.matD\tempCurve.vecG) * tempCurve.deltaNormVals );
		tempCurve.vecDDeltaVals = tempCurve.vecDeltaVals(:,2:end) - tempCurve.vecDeltaVals(:,1:end-1);
		tempCurve.alphaVals = (-(tempCurve.matD\tempCurve.vecG)' * tempCurve.vecDDeltaVals) ...
		 ./ ( norm(tempCurve.matD\tempCurve.vecG) * sqrt(sum(tempCurve.vecDDeltaVals.^2,1)) );
		%
		levMarqModCurve = tempCurve;
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
	  levCurve.sVals, levCurve.pVals, 'o-', ...
	  levMarqCurve.sVals, levMarqCurve.pVals, 'x-', ...
	  levMarqModCurve.sVals, levMarqModCurve.pVals, '^-', ...
	  s0, p0, "s", "linewidth", 3, "markersize", 20, "color", [0.9,0.5,0.5], ...
	  s0, p0, "+", "linewidth", 3, "markersize", 20, "color", [0.9,0.5,0.5],
	  s1, p1, "s", "linewidth", 3, "markersize", 20, "color", [0.5,0.5,1.0], ...
	  s1, p1, "x", "linewidth", 3, "markersize", 20, "color", [0.5,0.5,1.0] );
	hold off;
	%
	%
	%
	numFigs++; figure(numFigs);
	loglog( ...
	  1.0./levCurve.lambdaVals(2:end-1) - 1.0, levCurve.omegaVals(2:end-1), 'o-', ...
	  1.0./levMarqCurve.lambdaVals(2:end-1) - 1.0, levMarqCurve.omegaVals(2:end-1), 'x-', ...
	  1.0./levMarqModCurve.lambdaVals(2:end-1) - 1.0, levMarqModCurve.omegaVals(2:end-1), '^-', ...
	  1.0./levCurve.lambdaVals(2:end-1) - 1.0, levCurve.omegaModelVals(2:end-1), 'o-', ...
	  1.0./levMarqCurve.lambdaVals(2:end-1) - 1.0, levMarqCurve.omegaModelVals(2:end-1), 'x-', ...
	  1.0./levMarqModCurve.lambdaVals(2:end-1) - 1.0, levMarqModCurve.omegaModelVals(2:end-1), '^-' );
	xlabel( "mu" );
	ylabel( "omega" );
	title( "omega vs mu" );
	legend( ...
	  "Lev", ...
	  "LevMarq", ...
	  "LevMarqMod", ...
	  "Lev (model)", ...
	  "LevMarq (model)", ...
	  "LevMarqMod (model)", ...
	  "location", "SouthEast" );
	grid on;
	%
	if ( mygetfield(prm,"showMuLimits",false) )
		numFigs++; figure(numFigs);
		loglog( ...
		  1.0./(eps+levCurve.lambdaVals) - 1.0 + 2*eps, levCurve.omegaVals, 'o-', ...
		  1.0./(eps+levMarqCurve.lambdaVals) - 1.0 + 2*eps, levMarqCurve.omegaVals, 'x-',...
		  1.0./(eps+levMarqModCurve.lambdaVals) - 1.0 + 2*eps, levMarqModCurve.omegaVals, '^-',...
		  1.0./(eps+levCurve.lambdaVals) - 1.0 + 2*eps, levCurve.omegaModelVals, 'o-', ...
		  1.0./(eps+levMarqCurve.lambdaVals) - 1.0 + 2*eps, levMarqCurve.omegaModelVals, 'x-', ...
		  1.0./(eps+levMarqModCurve.lambdaVals) - 1.0 + 2*eps, levMarqModCurve.omegaModelVals, '^-' );
		xlabel( "mu" );
		ylabel( "omega" );
		title( "omega vs mu" );
		legend( ...
	  "Lev", ...
	  "LevMarq", ...
	  "LevMarqMod", ...
	  "Lev (model)", ...
	  "LevMarq (model)", ...
	  "LevMarqMod (model)", ...
		  "location", "NorthEast" );
		grid on;
	end
	%
	%
	%
	numFigs++; figure(numFigs);
	semilogx( ...
	  levCurve.muVals(2:end-1), levCurve.deltaNormVals(2:end-1), 'o-', ...
	  levMarqCurve.muVals(2:end-1), levMarqCurve.deltaNormVals(2:end-1), 'x-', ...
	  levMarqModCurve.muVals(2:end-1), levMarqModCurve.deltaNormVals(2:end-1), '^-' );
	xlabel( "mu" );
	ylabel( "||delta||" );
	title( "||delta|| vs mu" );
	legend( ...
	  "Lev", ...
	  "LevMarq", ...
	  "LevMarqMod", ...
	  "location", "SouthEast" );
	grid on;
	%
	%
	%
	return;
	
	
	
	%
	numFigs++; figure(numFigs);
	semilogx( ...
	  levCurve.muVals(2:end), levCurve.alphaVals, 'o-', ...
	  levMarqCurve.muVals(2:end), levMarqCurve.alphaVals, 'x-' );
	xlabel( "mu" );
	ylabel( "alpha" );
	title( "alpha vs mu" );
	legend( ...
	  "Lev", ...
	  "LevMarq", ...
	  "location", "SouthEast" );
	grid on;
	%
	%
	%
	numFigs++; figure(numFigs);
	semilogx( ...
	  cent(levCurve.muVals(2:end)), diff(levCurve.alphaVals)./diff(levCurve.muVals(2:end)), 'o-', ...
	  cent(levMarqCurve.muVals(2:end)), diff(levMarqCurve.alphaVals)./diff(levMarqCurve.muVals(2:end)), 'x-' );
	xlabel( "mu" );
	ylabel( "d alpha / d mu" );
	title( "d alpha / d mu vs mu" );
	legend( ...
	  "Lev", ...
	  "LevMarq", ...
	  "location", "SouthEast" );
	grid on;
	%
	return;
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( ...
	  levCurve.lambdaVals, levCurve.omegaVals, 'o-', ...
	  levMarqCurve.lambdaVals, levMarqCurve.omegaVals, 'x-' );
	xlabel( "lambda" );
	ylabel( "omega" );
	title( "omega vs lambda" );
	legend( ...
	  "Lev", ...
	  "LevMarq", ...
	  "location", "NorthEast" );
	grid on;
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( ...
	  levCurve.lambdaVals, levCurve.deltaNormVals, 'o-', ...
	  levMarqCurve.lambdaVals, levMarqCurve.deltaNormVals, 'x-' );
	xlabel( "lambda" );
	ylabel( "||delta||" );
	title( "||delta|| vs lambda" );
	legend( ...
	  "Lev", ...
	  "LevMarq", ...
	  "location", "NorthWest" );
	grid on;
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( ...
	  levCurve.lambdaVals, levCurve.dacVals, 'o-', ...
	  levMarqCurve.lambdaVals, levMarqCurve.dacVals, 'x-' );
	xlabel( "lambda" );
	ylabel( "DAC" );
	title( "DAC vs lambda" );
	legend( ...
	  "Lev", ...
	  "LevMarq", ...
	  "location", "NorthWest" );
	grid on;
return;
end
