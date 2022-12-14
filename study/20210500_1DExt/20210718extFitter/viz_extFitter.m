	clear;
	commondefs;
	thisFile = "viz_extFitter";
	numFigs = 0;
	setprngstates(0);
	%
	secret_p = 3.23
	secret_xExt = 0.14
	secret_fExt = 0.024
	secret_f1 = 0.76;
	funchF = @(x)( secret_fExt + secret_f1*abs(x-secret_xExt).^secret_p );
	%%%funchF = @(x)( (x<0).*abs(x).^4.0 + (x>0).*abs(x).^0.5 );
	%
	numPts = 20;
	%xVals = 0.3+sort(abs(randn(1,numPts)));
	xVals = sort(randn(1,numPts));
	%
	fVals = funchF(xVals);
	%fVals .*= 1.0 + 0.01*randn(1,numPts);
	%
	%
	[ foo, nC ] = min(abs(fVals));
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
