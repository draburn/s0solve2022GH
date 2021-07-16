	clear;
	commondefs;
	thisFile = "test_extFitterOneSided";
	numFigs = 0;
	%
	%
	setprngstates(0);
	xVals = sort(abs(randn(1,10)));
	fVals = abs(xVals).^4;
	%
	%
	nOfClosestPt = 1;
	%
	%
	%sVals = linspace(-0.1,0.0,31);
	%pVals = linspace(3.8,4.2,51);
	sVals = linspace(-0.2,0.2,51);
	pVals = linspace(3.0,5.0,51);
	[ sMesh, pMesh ] = meshgrid( sVals, pVals );
	size1 = size(sMesh,1);
	size2 = size(sMesh,2);
	[ omegaMesh, bigG0Mesh, bigG1Mesh, dat ] = extFitterOneSided_getOmegaMesh( ...
	  sMesh, pMesh, xVals, fVals, nOfClosestPt );
	xExtMesh = dat.bigX0 + dat.bigX1*sMesh;
	fExtMesh = dat.bigF0 + dat.bigF1*bigG0Mesh;
	sqrtOmegaMesh = sqrt(omegaMesh);
	%
	%
	%
	meshZ = sqrtOmegaMesh;
	strZ = "sqrt(omega)";
	numFigs++; figure(numFigs);
	plot( sMesh((size1+[-3,1,3])/2,:)', meshZ((size1+[-3,1,3])/2,:)', 'o-' );
	grid on;
	xlabel( "s" );
	ylabel( strZ );
	title( [ strZ " vs s" ] );
	%
	numFigs++; figure(numFigs);
	plot( pMesh(:,(size2+[-3,1,3])/2), meshZ(:,(size2+[-3,1,3])/2), 'o-' );
	grid on;
	xlabel( "p" );
	ylabel( strZ );
	title( [ strZ " vs p" ] );
	%
	numContours = 30;
	numColors = numContours+1;
	numFigs++; figure(numFigs);
	contourf( sVals, pVals, meshZ, numContours );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( [ strZ " vs s, p" ] );
	%
	%
	%
	meshZ = bigG0Mesh;
	strZ = "G0";
	numFigs++; figure(numFigs);
	plot( sMesh((size1+[-3,1,3])/2,:)', meshZ((size1+[-3,1,3])/2,:)', 'o-' );
	grid on;
	xlabel( "s" );
	ylabel( strZ );
	title( [ strZ " vs s" ] );
	%
	numFigs++; figure(numFigs);
	plot( pMesh(:,(size2+[-3,1,3])/2), meshZ(:,(size2+[-3,1,3])/2), 'o-' );
	grid on;
	xlabel( "p" );
	ylabel( strZ );
	title( [ strZ " vs p" ] );
	%
	numContours = 30;
	numColors = numContours+1;
	numFigs++; figure(numFigs);
	contourf( sVals, pVals, meshZ, numContours );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( [ strZ " vs s, p" ] );
	%
	%
	%
	meshZ = bigG1Mesh;
	strZ = "G1";
	numFigs++; figure(numFigs);
	plot( sMesh((size1+[-3,1,3])/2,:)', meshZ((size1+[-3,1,3])/2,:)', 'o-' );
	grid on;
	xlabel( "s" );
	ylabel( strZ );
	title( [ strZ " vs s" ] );
	%
	numFigs++; figure(numFigs);
	plot( pMesh(:,(size2+[-3,1,3])/2), meshZ(:,(size2+[-3,1,3])/2), 'o-' );
	grid on;
	xlabel( "p" );
	ylabel( strZ );
	title( [ strZ " vs p" ] );
	%
	numContours = 30;
	numColors = numContours+1;
	numFigs++; figure(numFigs);
	contourf( sVals, pVals, meshZ, numContours );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( [ strZ " vs s, p" ] );
	%
	%
	%
	meshZ = xExtMesh;
	strZ = "x_{ext}";
	numFigs++; figure(numFigs);
	plot( sMesh((size1+[-3,1,3])/2,:)', meshZ((size1+[-3,1,3])/2,:)', 'o-' );
	grid on;
	xlabel( "s" );
	ylabel( strZ );
	title( [ strZ " vs s" ] );
	%
	numFigs++; figure(numFigs);
	plot( pMesh(:,(size2+[-3,1,3])/2), meshZ(:,(size2+[-3,1,3])/2), 'o-' );
	grid on;
	xlabel( "p" );
	ylabel( strZ );
	title( [ strZ " vs p" ] );
	%
	numContours = 30;
	numColors = numContours+1;
	numFigs++; figure(numFigs);
	contourf( sVals, pVals, meshZ, numContours );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( [ strZ " vs s, p" ] );
	%
	%
	%
	meshZ = fExtMesh;
	strZ = "f_{ext}";
	numFigs++; figure(numFigs);
	plot( sMesh((size1+[-3,1,3])/2,:)', meshZ((size1+[-3,1,3])/2,:)', 'o-' );
	grid on;
	xlabel( "s" );
	ylabel( strZ );
	title( [ strZ " vs s" ] );
	%
	numFigs++; figure(numFigs);
	plot( pMesh(:,(size2+[-3,1,3])/2), meshZ(:,(size2+[-3,1,3])/2), 'o-' );
	grid on;
	xlabel( "p" );
	ylabel( strZ );
	title( [ strZ " vs p" ] );
	%
	numContours = 30;
	numColors = numContours+1;
	numFigs++; figure(numFigs);
	contourf( sVals, pVals, meshZ, numContours );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "s" );
	ylabel( "p" );
	title( [ strZ " vs s, p" ] );
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( xVals, fVals, 'ko-', 'linewidth', 5, 'markersize', 20 );
	grid on;
	xlabel( "x" );
	ylabel( "f" );
	hold on;
	yVals = ( xVals - dat.bigX0 ) / dat.bigX1;
	for i1=(size1+[-3,1,3])/2
	for i2=(size2+[-3,1,3])/2
		gModelVals = bigG0Mesh(i1,i2) ...
		  + bigG1Mesh(i1,i2)*(abs(yVals-sMesh(i1,i2)).^pMesh(i1,i2));
		fModelVals = dat.bigF0 + dat.bigF1*gModelVals;
		plot( xVals, fModelVals, 'x-', 'linewidth', 3, 'markersize', 10 );
	end
	end
	hold off;
	%
	%
	%
	numFigs++; figure(numFigs);
	semilogy( xVals, fVals, 'ko-', 'linewidth', 5, 'markersize', 20 );
	grid on;
	xlabel( "x" );
	ylabel( "f" );
	hold on;
	yVals = ( xVals - dat.bigX0 ) / dat.bigX1;
	for i1=(size1+[-3,1,3])/2
	for i2=(size2+[-3,1,3])/2
		gModelVals = bigG0Mesh(i1,i2) ...
		  + bigG1Mesh(i1,i2)*(abs(yVals-sMesh(i1,i2)).^pMesh(i1,i2));
		fModelVals = dat.bigF0 + dat.bigF1*gModelVals;
		semilogy( xVals, fModelVals, 'x-', 'linewidth', 3, 'markersize', 10 );
	end
	end
	hold off;
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( xVals, fVals-fVals, 'k-', 'linewidth', 2 );
	grid on;
	xlabel( "x" );
	ylabel( "f_{model} - f" );
	hold on;
	yVals = ( xVals - dat.bigX0 ) / dat.bigX1;
	for i1=(size1+[-3,1,3])/2
	for i2=(size2+[-3,1,3])/2
		gModelVals = bigG0Mesh(i1,i2) ...
		  + bigG1Mesh(i1,i2)*(abs(yVals-sMesh(i1,i2)).^pMesh(i1,i2));
		fModelVals = dat.bigF0 + dat.bigF1*gModelVals;
		plot( xVals, fModelVals-fVals, 'x-', 'linewidth', 3, 'markersize', 10 );
	end
	end
	hold off;
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( xVals, fVals./abs(fVals)-1.0, 'k-', 'linewidth', 2 );
	grid on;
	xlabel( "x" );
	ylabel( "f_{model}/|f| - 1" );
	hold on;
	yVals = ( xVals - dat.bigX0 ) / dat.bigX1;
	for i1=(size1+[-3,1,3])/2
	for i2=(size2+[-3,1,3])/2
		gModelVals = bigG0Mesh(i1,i2) ...
		  + bigG1Mesh(i1,i2)*(abs(yVals-sMesh(i1,i2)).^pMesh(i1,i2));
		fModelVals = dat.bigF0 + dat.bigF1*gModelVals;
		plot( xVals, fModelVals./abs(fVals)-1.0, 'x-', 'linewidth', 3, 'markersize', 10 );
	end
	end
	hold off;
return;
