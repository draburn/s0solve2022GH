	clear;
	commondefs;
	thisFile = "test_extFitterOneSided2";
	numFigs = 0;
	%
	%
	setprngstates(0);
	%xVals = 0.3+sort(abs(randn(1,20)));
	xVals = 1.0+sort(abs(randn(1,20)));
	%xVals = linspace(0.2,1.0,5);
	secret_xExt = 0.0;
	secret_fExt = 0.0;
	secret_p = 3.23;
	secret_f1 = 0.3;
	funchF = @(x)( secret_fExt + secret_f1*abs(x-secret_xExt).^secret_p );
	%
	
	fVals = funchF(xVals);
	nOfClosestPt = 1;
	
	doAuxPlots = true;
	if (doAuxPlots)
		sValsA = linspace(-4.2,0.1,5001);
		%sValsA = linspace(-0.76,-0.75,5001);
		pValsA = 3.23+(sValsA+0.75406)*(-1.85);
		pValsA2 = 3.23+(sValsA+0.75406)/(1.85);
		numFigs++; figure(numFigs);
		plot( sValsA, pValsA, 'o-', sValsA, pValsA2, 'x-' );
		grid on;
		xlabel( "s" );
		ylabel( "p" );
		axis( "equal" );
		[ omegaValsA, bigF0ValsA, bigF1ValsA, dat ] = extFitterOneSided_getOmegaMesh2( ...
		  sValsA, pValsA, xVals, fVals, nOfClosestPt );
		[ omegaValsA2, bigF0ValsA2, bigF1ValsA2, dat ] = extFitterOneSided_getOmegaMesh2( ...
		  sValsA, pValsA2, xVals, fVals, nOfClosestPt );
		numFigs++; figure(numFigs);
		plot( sValsA, sqrt(omegaValsA), 'o-', sValsA, sqrt(omegaValsA2), 'x-' );
		grid on;
		xlabel( "s" );
		ylabel( "sqrt(omega)" );
		
		sqrt(omegaValsA(1))/sqrt(omegaValsA2(1))
		
		%return;
	end
	
	
	%
	%
	numContours = 30;
	numColors = numContours+1;
	sVals = linspace(-5.0,0.0,51);
	pVals = linspace(0.1,10.0,51);

	
	sVals = linspace(-5.0,-2.0,51);
	pVals = linspace(4.0,10.0,51);
	
	if (0)
	sVals = linspace(-1.0,0.0,51);
	pVals = linspace(2.0,5.0,51);
	
	sVals = linspace(-0.2,-0.0,51);
	pVals = linspace(3.0,3.5,51);
	
	sVals = linspace(-0.12,-0.08,51);
	pVals = linspace(3.18,3.26,51);
	
	sVals = linspace(-0.106,-0.100,51);
	pVals = linspace(3.22,3.24,51);
	
	sVals = linspace(-0.1034,-0.1020,51);
	pVals = linspace(3.228,3.234,51);
	end
	
	if (0)
	sVals = linspace(-0.5,-3.0,51);
	pVals = linspace(2.0,6.0,51);
	
	sVals = linspace(-1.0,-0.5,51);
	pVals = linspace(2.5,4.0,51);
	
	sVals = linspace(-0.9,-0.6,51);
	pVals = linspace(3.0,3.6,51);
	% Converging!
	
	sVals = linspace(-0.78,-0.74,51);
	pVals = linspace(3.20,3.28,51);
	
	sVals = linspace(-0.76,-0.74,51);
	pVals = linspace(3.21,3.24,51);
	% Go tight...
	
	sVals = linspace(-0.752,-0.748,51);
	pVals = linspace(3.226,3.232,51);
	
	sVals = linspace(-0.7506,-0.7504,51);
	pVals = linspace(3.229,3.231,51);
	end
	
	%sVals = linspace(-1.0,-0.1,51);
	%pVals = linspace(2.0,6.0,51);
	
	%sVals = linspace(-0.5,-0.1,51);
	%pVals = linspace(2.5,4.0,51);
	
	%sVals = linspace(-0.35,-0.25,51);
	%pVals = linspace(3.1,3.4,51);
	
	%sVals = linspace(-0.75047,-0.75045,51);
	%pVals = linspace(3.22999,3.23001,51);
	
	%
	%sVals = linspace(-1.0,-0.3,51);
	%pVals = linspace(2.0,4.0,51);
	%pVals = linspace(2.0,5.0,51);
	%sVals = linspace(-1.0,-0.4,51);
	%sVals = linspace(-0.7505,-0.7504,101);
	%pVals = linspace(5.2299,5.2301,101);
	[ sMesh, pMesh ] = meshgrid( sVals, pVals );
	size1 = size(sMesh,1);
	size2 = size(sMesh,2);
	[ omegaMesh, bigF0Mesh, bigF1Mesh, dat ] = extFitterOneSided_getOmegaMesh2( ...
	  sMesh, pMesh, xVals, fVals, nOfClosestPt );
	xExtMesh = dat.bigX0 + dat.bigX1*sMesh;
	[ omegaCrit, i1Crit, i2Crit ] = minmin(omegaMesh);
	if (doAuxPlots)
		echo__omegaCrit = omegaCrit
		i1Crit = 40
		i2Crit = 19
		omegaCrit = omegaMesh(i1Crit,i2Crit)
	end
	sCrit = sMesh(i1Crit,i2Crit)
	pCrit = pMesh(i1Crit,i2Crit)
	bigF0Crit = bigF0Mesh(i1Crit,i2Crit)
	bigF1Crit = bigF1Mesh(i1Crit,i2Crit)
	xExtCrit = xExtMesh(i1Crit,i2Crit)
	%
	%
	%
	%meshZ = omegaMesh.^0.25; strZ = "omega^{0.25}";
	meshZ = omegaMesh.^0.5; strZ = "omega^{0.5}";
	numFigs++; figure(numFigs);
	semilogy( sMesh(i1Crit+[-1,0,1],:)', meshZ(i1Crit+[-1,0,1],:)', 'o-' );
	grid on;
	xlabel( "s" );
	ylabel( strZ );
	title( [ strZ " vs s" ] );
	%
	numFigs++; figure(numFigs);
	semilogy( pMesh(:,i2Crit+[-1,0,1]), meshZ(:,i2Crit+[-1,0,1]), 'o-' );
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
	if (doAuxPlots)
		plot( sValsA, pValsA, 'o-', sValsA, pValsA2, 'x-' );
	end
	plot( sCrit, pCrit, 'wo', 'linewidth', 3, 'markersize', 25 );
	hold off;
	%
	%
	%
	xLo = min([ xExtCrit, min(xVals) ]);
	xLo = max([ xLo, secret_xExt - 0.1 ]);
	xHi = max([ xExtCrit, max(xVals) ]);
	viz_xVals = linspace(xLo-0.1*(xHi-xLo),xHi+0.1*(xHi-xLo),1001);
	viz_fVals = funchF(viz_xVals);
	funchG = @(x)( bigF0Crit + bigF1Crit*( abs((x-dat.bigX0)/dat.bigX1-sCrit).^pCrit ) );
	gVals = funchG(xVals);
	viz_gVals = funchG(viz_xVals);
	%
	numFigs++; figure(numFigs);
	plot( ...
	  viz_xVals, viz_gVals, 'x-', 'linewidth', 1, 'markersize', 10, ...
	  viz_xVals, viz_fVals, '+-', 'linewidth', 1, 'markersize', 10, ...
	  xExtCrit, bigF0Crit, 'x', 'linewidth', 3, 'markersize', 35, ...
	  secret_xExt, secret_fExt, '+', 'linewidth', 3, 'markersize', 35, ...
	  xVals, gVals, 'x', 'linewidth', 4, 'markersize', 20, ...
	  xVals, fVals, 'k+', 'linewidth', 4, 'markersize', 20 );
	grid on;
	xlabel( "x" );
	ylabel( "f" );
	title( "f vs x" );
	%
	numFigs++; figure(numFigs);
	plot( ...
	  viz_xVals, viz_gVals-viz_fVals, 'x-', 'linewidth', 1, 'markersize', 10, ...
	  xVals, gVals-fVals, 'x', 'linewidth', 4, 'markersize', 20 );
	grid on;
	xlabel( "x" );
	ylabel( "g - f" );
	title( "g - f vs x" );
	%
	numFigs++; figure(numFigs);
	plot( ...
	  viz_xVals, viz_gVals, 'x-', 'linewidth', 1, 'markersize', 10, ...
	  xVals, fVals, 'k+', 'linewidth', 5, 'markersize', 30 );
	grid on;
	xlabel( "x" );
	ylabel( "f" );
	title( "f vs x" );
	%
	%
	
	
	%
	%
	%
	meshZ = (pMesh-1.0)./(2.0*sMesh); strZ = "alpha thing";
	meshZ = abs(meshZ+1.04);
	numFigs++; figure(numFigs);
	semilogy( sMesh(i1Crit+[-1,0,1],:)', meshZ(i1Crit+[-1,0,1],:)', 'o-' );
	grid on;
	xlabel( "s" );
	ylabel( strZ );
	title( [ strZ " vs s" ] );
	%
	numFigs++; figure(numFigs);
	semilogy( pMesh(:,i2Crit+[-1,0,1]), meshZ(:,i2Crit+[-1,0,1]), 'o-' );
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
	if (doAuxPlots)
		plot( sValsA, pValsA, 'o-', sValsA, pValsA2, 'x-' );
	end
	plot( sCrit, pCrit, 'wo', 'linewidth', 3, 'markersize', 25 );
	hold off;

