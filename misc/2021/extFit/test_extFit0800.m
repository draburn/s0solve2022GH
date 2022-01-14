	clear;
	commondefs;
	thisFile = "test_extFit0800";
	numFigs = 0;
	setprngstates(0);
	%
	%secret_p = 3.23
	%secret_xExt = 0.14
	%secret_fExt = 0.024
	%secret_f1 = 0.76;
	secret_xExt = 0.0
	secret_p = 6.73
	secret_fExt = 0.0
	secret_f1 = 1.0;
	funchF = @(x)( secret_fExt + secret_f1*abs(x-secret_xExt).^secret_p );
	%%%funchF = @(x)( (x<0).*abs(x).^4.0 + (x>0).*abs(x).^0.5 );
	%
	numPts = 10;
	%xVals = randn(1,numPts);
	xVals = 1.0+abs(randn(1,numPts));
	%
	xVals = sort(xVals);
	numPts = max(size(xVals));
	fVals = funchF(xVals);
	fVals .*= 1.0 + 0.01*randn(1,numPts);
	%
	%
	%wVals = ones(size(fVals));
	wVals = 1.0./(fVals.^2+eps*min(fVals.^2));
	wVals /= sum(wVals);
	tic
	prm_viz = [];
	%prm.sLo = -0.2;
	%prm.sHi = 0.2;
	%prm.pLo = 2.0;
	%prm.pHi = 6.0;
	s0 = 0.0; p0 = 3.0;
	%
	prm_viz.numFigs0 = numFigs+10;
	s = s0; p = p0;
	viz_extFitPt( xVals, fVals, s, p, wVals, prm_viz );
	%
	prm_findFit = [];
	prm_findFit.sMin = min(xVals) - 1.0*(max(xVals)-min(xVals));
	prm_findFit.sMax = max(xVals) + 1.0*(max(xVals)-min(xVals));
	prm_findFit.prm_findStep.useLevMarq = true;
	[ s1, p1, retCode ] = extFit__findFit( ...
	  s0, p0, xVals, fVals, wVals, prm_findFit );
	msg( thisFile, __LINE__, sprintf( ...
	  "extFit__findFit() returned %s.", retcode2str(retCode) ) );
	%
	prm_calcAboutPt = [];
	[ rhoVals0, bigF00, bigF10, omega0, vecG0, matH0 ] = extFit__calcAboutPt( ...
	  s0, p0, xVals, fVals, wVals, prm_calcAboutPt );
	[ rhoVals1, bigF01, bigF11, omega1, vecG1, matH1 ] = extFit__calcAboutPt( ...
	  s1, p1, xVals, fVals, wVals, prm_calcAboutPt );
	funchFModel0 = @(x)(  bigF00 + bigF10 * abs( x - s0 ).^p0  );
	funchFModel1 = @(x)(  bigF01 + bigF11 * abs( x - s1 ).^p1  );
	%
	viz_numPts = 1000;
	viz_xLo = min(xVals) - 0.3 * ( max(xVals) - min(xVals) );
	viz_xHi = max(xVals) + 0.3 * ( max(xVals) - min(xVals) );
	viz_xVals = linspace( viz_xLo, viz_xHi, viz_numPts );
	%
	numFigs++; figure(numFigs);
	plot( ...
	  xVals, funchFModel1(xVals), 'x', 'linewidth', 3, 'markersize', 15, ...
	  viz_xVals, funchFModel1(viz_xVals), 'x-', 'linewidth', 1, 'markersize', 5, ...
	  xVals, funchFModel0(xVals), 'o', 'linewidth', 3, 'markersize', 15, ...
	  viz_xVals, funchFModel0(viz_xVals), 'o-', 'linewidth', 1, 'markersize', 5, ...
	  xVals, funchF(xVals), '+', 'linewidth', 3, 'markersize', 15, ...
	  viz_xVals, funchF(viz_xVals), '+-', 'linewidth', 1, 'markersize', 5, ...
	  xVals, fVals, 'ko', 'linewidth', 5, 'markersize', 25, ...
	  xVals, fVals, 'k+', 'linewidth', 5, 'markersize', 25 );
	grid on;
	xlabel( "X" );
	ylabel( "f" );
	title( "f vs x" );
	legend( ...
	  "model 1 pts", ...
	  "model 1", ...
	  "model 0 pts", ...
	  "model 0", ...
	  "base pts", ...
	  "base", ...
	  "actual pts", ...
	  "actual pts", ...
	  "location", "north" );
	%
	numFigs++; figure(numFigs);
	plot( ...
	  xVals, funchFModel1(xVals)-funchF(xVals), 'x', 'linewidth', 3, 'markersize', 15, ...
	  viz_xVals, funchFModel1(viz_xVals)-funchF(viz_xVals), 'x-', 'linewidth', 1, 'markersize', 5, ...
	  xVals, fVals-funchF(xVals), 'ko', 'linewidth', 5, 'markersize', 25, ...
	  xVals, fVals-funchF(xVals), 'k+', 'linewidth', 5, 'markersize', 25 );
	grid on;
	xlabel( "X" );
	ylabel( "f - f base" );
	title( "f - f base vs x" );
	legend( ...
	  "model 1 pts", ...
	  "model 1", ...
	  "actual pts", ...
	  "actual pts", ...
	  "location", "north" );
	%
	numFigs++; figure(numFigs);
	plot( ...
	  xVals, funchFModel1(xVals)-fVals, 'x-', 'linewidth', 3, 'markersize', 15, ...
	  xVals, fVals-fVals, 'ko-', 'linewidth', 5, 'markersize', 25, ...
	  xVals, fVals-fVals, 'k+-', 'linewidth', 5, 'markersize', 25 );
	grid on;
	xlabel( "X" );
	ylabel( "f - f pts" );
	title( "f - f pts vs x" );
	legend( ...
	  "model 1 pts", ...
	  "actual pts", ...
	  "actual pts", ...
	  "location", "north" );
	%
	return
	
	
	for n=1:15
		[ s, p ] = extFit__findStep( s, p, xVals, fVals, nC, [] );
	end
	toc
	return;
	
	
	nFit = nC;
	if ( 1==nC )
		nL = 1;
		nR = 2;
		sHi = xVals(2);
		sLo = 2.0*xVals(1)-xVals(numPts);
	elseif (numPts==nC)
		nL = numPts-1;
		nR = numPts;
		sLo = xVals(numPts-1);
		sHi = 2.0*xVals(numPts)-xVals(1);
	else
		sLo = xVals(nC-1);
		sHi = xVals(nC+1);
		if (0)
		% DRaburn 2021.07.18.
		% No, we can't use quad fit;
		% For twoPt fit, we MUST have points on either side of s.
		% We could do this by adjusting the points based on s,
		% but, for fixed s, we must step out to either side.
		vecX = xVals(nC-1:nC+1)';
		vecF = fVals(nC-1:nC+1)';
		matX = [ ones(3,1), vecX, vecX.^2 ];
		vecC = matX\vecF;
		xQuad = -vecC(2)/(2.0*vecC(3));
		if ( xQuad > xVals(nC) )
			nL = nC;
			nR = nC+1;
		else
			nL = nC-1;
			nR = nC;
		end
		end
		nL = nC-1;
		nR = nC+1;
	end
	msg( thisFile, __LINE__, sprintf( ...
	  "nFit = %d, nL = %d, nR = %d, numPts = %d.", nFit, nL, nR, numPts ) );
	%
	numContours = 30;
	numColors = numContours+1;
	sMin = sLo-0.01*(sHi-sLo);
	sMax = sHi+0.01*(sHi-sLo);
	sVals = linspace(sMin,sMax,51);
	pVals = linspace(0.1,10.0,50);
	[ sMesh, pMesh ] = meshgrid( sVals, pVals );
	%
	[ bigF0Mesh, bigF1Mesh, omegaMesh ] = extFit__calcMesh( ...
	  sMesh, pMesh, xVals, fVals, nC );
	meshZ = omegaMesh.^0.5; strZ = "omega^{0.5}";
	numFigs++; figure(numFigs);
	contourf( sVals, pVals, meshZ, numContours );
	axis("square");
	%axis("equal");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( [ strZ " vs s, p" ] );
	
	return;
	%
	
	vecX = [ 0.0; 3.0 ];
	prm.calcH2 = true;
	for n=1:10
	[ bigF0, bigF1, rhoVals, omega, vecG, matH, matH2 ] = extFit__calcAboutPt( ...
	  vecX(1), vecX(2), xVals, fVals, nC, [], prm )
	vecX -= matH\vecG
	end
	return;
	
	
	nFit = nC;
	if ( 1==nC )
		nL = 1;
		nR = 2;
		sHi = xVals(2);
		sLo = 2.0*xVals(1)-xVals(numPts);
	elseif (numPts==nC)
		nL = numPts-1;
		nR = numPts;
		sLo = xVals(numPts-1);
		sHi = 2.0*xVals(numPts)-xVals(1);
	else
		sLo = xVals(nC-1);
		sHi = xVals(nC+1);
		if (0)
		% DRaburn 2021.07.18.
		% No, we can't use quad fit;
		% For twoPt fit, we MUST have points on either side of s.
		% We could do this by adjusting the points based on s,
		% but, for fixed s, we must step out to either side.
		vecX = xVals(nC-1:nC+1)';
		vecF = fVals(nC-1:nC+1)';
		matX = [ ones(3,1), vecX, vecX.^2 ];
		vecC = matX\vecF;
		xQuad = -vecC(2)/(2.0*vecC(3));
		if ( xQuad > xVals(nC) )
			nL = nC;
			nR = nC+1;
		else
			nL = nC-1;
			nR = nC;
		end
		end
		nL = nC-1;
		nR = nC+1;
	end
	msg( thisFile, __LINE__, sprintf( ...
	  "nFit = %d, nL = %d, nR = %d, numPts = %d.", nFit, nL, nR, numPts ) );
	%
	numContours = 30;
	numColors = numContours+1;
	sMin = sLo-0.01*(sHi-sLo);
	sMax = sHi+0.01*(sHi-sLo);
	sVals = linspace(sMin,sMax,50);
	pVals = linspace(0.1,10.0,50);
	%sVals = linspace(xVals(nC)-0.1,xVals(nC)+0.1,101);
	%pVals = linspace(4.0,6.0,21);
	%
	[ sMesh, pMesh ] = meshgrid( sVals, pVals );
	size1 = size(sMesh,1);
	size2 = size(sMesh,2);
	%
	if (1)
		%
		[ sSolved, pSolved ] = extFit_iterSolve( 0.1, 2.0, xVals, fVals, nFit )
		%
		for i1=1:size1
		for i2=1:size2
			[ omega, bigF0, bigF1 ] = extFit_getOmega( ...
			  sMesh(i1,i2), pMesh(i1,i2), xVals, fVals, nFit );
			omegaMesh1(i1,i2) = omega;
			bigF0Mesh1(i1,i2) = bigF0;
			bigF1Mesh1(i1,i2) = bigF1;
		end
		end
		dat1.bigX0 = 0.0;
		dat1.bigX1 = 1.0;
	else
	[ omegaMesh1, bigF0Mesh1, bigF1Mesh1, dat1 ] = extFitter_onePtFit( ...
	  sMesh, pMesh, xVals, fVals, nFit );
	end
	xExtMesh1 = dat1.bigX0 + dat1.bigX1*sMesh;
	[ omegaCrit1, i1Crit1, i2Crit1 ] = minmin(omegaMesh1);
	sCrit1 = sMesh(i1Crit1,i2Crit1);
	pCrit1 = pMesh(i1Crit1,i2Crit1)
	xExtCrit1 = xExtMesh1(i1Crit1,i2Crit1)
	fExtCrit1 = bigF0Mesh1(i1Crit1,i2Crit1)
	bigF0Crit1 = bigF0Mesh1(i1Crit1,i2Crit1);
	bigF1Crit1 = bigF1Mesh1(i1Crit1,i2Crit1);
	pLMesh = pMesh;
	pRMesh = pMesh;
	[ omegaMesh2, bigF0Mesh2, bigFLMesh2, bigFRMesh2, dat2 ] = extFitter_twoPtFit( ...
	  sMesh, pLMesh, pRMesh, xVals, fVals, nL, nR );
	xExtMesh2 = dat1.bigX0 + dat1.bigX1*sMesh;
	[ omegaCrit2, i1Crit2, i2Crit2 ] = minmin(omegaMesh2);
	sCrit2 = sMesh(i1Crit2,i2Crit2);
	pLCrit2 = pLMesh(i1Crit2,i2Crit2)
	pRCrit2 = pRMesh(i1Crit2,i2Crit2)
	xExtCrit2 = xExtMesh2(i1Crit2,i2Crit2)
	fExtCrit2 = bigF0Mesh2(i1Crit2,i2Crit2)
	bigF0Crit2 = bigF0Mesh2(i1Crit2,i2Crit2);
	bigFLCrit2 = bigFLMesh2(i1Crit2,i2Crit2);
	bigFRCrit2 = bigFRMesh2(i1Crit2,i2Crit2);
	%
	%
	%
	meshZ = omegaMesh1.^0.5; strZ = "omega1^{0.5}";
	%meshZ = log(omegaMesh1); strZ = "log(omega1)";
	i1CritVals1 = cap( i1Crit1+[-1,0,1], 1, size1 );
	numFigs++; figure(numFigs);
	plot( sMesh(i1CritVals1,:)', meshZ(i1CritVals1,:)', 'o-' );
	%semilogy( sMesh(i1Crit1+[-1,0,1],:)', meshZ(i1Crit1+[-1,0,1],:)', 'o-' );
	grid on;
	xlabel( "s" );
	ylabel( strZ );
	title( [ strZ " vs s" ] );
	%
	numFigs++; figure(numFigs);
	i2CritVals1 = cap( i1Crit1+[-1,0,1], 1, size2 );
	plot( pMesh(:,i2CritVals1), meshZ(:,i2CritVals1), 'o-' );
	%semilogy( pMesh(:,i2Crit1+[-1,0,1]), meshZ(:,i2Crit1+[-1,0,1]), 'o-' );
	grid on;
	xlabel( "p" );
	ylabel( strZ );
	title( [ strZ " vs p" ] );
	%
	numFigs++; figure(numFigs);
	contourf( sVals, pVals, meshZ, numContours );
	axis("square");
	%axis("equal");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( [ strZ " vs s, p" ] );
	hold on;
	plot( sCrit1, pCrit1, 'wo', 'linewidth', 3, 'markersize', 25 );
	hold off;
	%
	%
	%
	meshZ = omegaMesh2.^0.5; strZ = "omega2^{0.5}";
	%meshZ = log(omegaMesh2); strZ = "log(omega2)";
	numFigs++; figure(numFigs);
	i1CritVals2 = cap( i1Crit2+[-1,0,1], 1, size1 );
	plot( sMesh(i1CritVals2,:)', meshZ(i1CritVals2,:)', 'o-' );
	%semilogy( sMesh(i1CritVals2,:)', meshZ(i1CritVals2,:)', 'o-' );
	grid on;
	xlabel( "s" );
	ylabel( strZ );
	title( [ strZ " vs s" ] );
	%
	numFigs++; figure(numFigs);
	i2CritVals2 = cap( i2Crit2+[-1,0,1], 1, size2 );
	plot( pMesh(:,i2CritVals2), meshZ(:,i2CritVals2), 'o-' );
	%semilogy( pMesh(:,i2CritVals2), meshZ(:,i2CritVals2), 'o-' );
	grid on;
	xlabel( "p" );
	ylabel( strZ );
	title( [ strZ " vs p" ] );
	%
	numFigs++; figure(numFigs);
	contourf( sVals, pVals, meshZ, numContours );
	axis("square");
	%axis("equal");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( [ strZ " vs s, p" ] );
	hold on;
	plot( sCrit2, (pLCrit2+pRCrit2)/2.0, 'wo', 'linewidth', 3, 'markersize', 25 );
	hold off;
	%
	%
	%
	xLo = min([ xExtCrit1, xExtCrit2, min(xVals) ]);
	xHi = max([ xExtCrit1, xExtCrit2, max(xVals) ]);
	viz_xVals = linspace(xLo-0.1*(xHi-xLo),xHi+0.1*(xHi-xLo),1001);
	viz_fVals = funchF(viz_xVals);
	%
	funchG1 = @(x)( bigF0Crit1 + bigF1Crit1*( abs((x-dat1.bigX0)/dat1.bigX1-sCrit1).^pCrit1 ) );
	gVals1 = funchG1(xVals);
	viz_gVals1 = funchG1(viz_xVals);
	%
	funchG2 = @(x)( bigF0Crit2 ...
	 + bigFLCrit2 * (x<xExtCrit2) .* abs((x-dat2.bigX0)/dat2.bigX1-sCrit2).^pLCrit2 ...
	 + bigFRCrit2 * (x>xExtCrit2) .* abs((x-dat2.bigX0)/dat2.bigX1-sCrit2).^pRCrit2 );
	gVals2 = funchG2(xVals);
	viz_gVals2 = funchG2(viz_xVals);
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( ...
	  viz_xVals, viz_gVals1, 'x-', 'linewidth', 1, 'markersize', 10, ...
	  viz_xVals, viz_gVals2, '+-', 'linewidth', 1, 'markersize', 10, ...
	  xExtCrit1, fExtCrit1, 'x', 'linewidth', 3, 'markersize', 35, ...
	  xExtCrit2, fExtCrit2, '+', 'linewidth', 3, 'markersize', 35, ...
	  secret_xExt, secret_fExt, 'ko', 'linewidth', 3, 'markersize', 35, ...
	  viz_xVals, viz_fVals, 'ko-', 'linewidth', 1, 'markersize', 10, ...
	  xVals, fVals, 'ko', 'linewidth', 4, 'markersize', 20 );
	grid on;
	xlabel( "x" );
	ylabel( "f" );
	title( "f vs x" );
	%
	numFigs++; figure(numFigs);
	plot( ...
	  viz_xVals, viz_gVals1-viz_fVals, 'x-', 'linewidth', 1, 'markersize', 10, ...
	  viz_xVals, viz_gVals2-viz_fVals, '+-', 'linewidth', 1, 'markersize', 10, ...
	  xVals, gVals1-fVals, 'x', 'linewidth', 4, 'markersize', 20, ...
	  xVals, gVals2-fVals, '+', 'linewidth', 4, 'markersize', 20 );
	grid on;
	xlabel( "x" );
	ylabel( "g - f" );
	title( "g - f vs x" );
