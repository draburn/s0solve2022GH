	clear;
	commondefs;
	thisFile = "test_extFitterOneSided";
	numFigs = 0;
	%
	%
	setprngstates(0);
	xVals = sort(abs(randn(1,200)));
	fVals = abs(xVals).^4;
	%
	%
	nOfClosestPt = 1;
	%
	%
	%sVals = linspace(-0.1,0.0,31);
	%pVals = linspace(3.8,4.2,51);
	sVals = linspace(-0.2,0.1,51);
	pVals = linspace(3.0,5.0,51);
	[ sMesh, pMesh ] = meshgrid( sVals, pVals );
	size1 = size(sMesh,1);
	size2 = size(sMesh,2);
	[ omegaMesh, bigG0Mesh, bigG1Mesh, dat ] = extFitterOneSided_getOmegaMesh( ...
	  sMesh, pMesh, xVals, fVals, nOfClosestPt );
	sqrtOmegaMesh = sqrt(omegaMesh);
	%
	numFigs++; figure(numFigs);
	plot( xVals, fVals, 'ko-', 'linewidth', 5, 'markersize', 20 );
	grid on;
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
return;

	  



function extFinder_viz( xVals, fVals, bigSLo, bigSHi, bigPLo, bigPHi, prm=[] )
	
	commondefs;
	thisFile = "extFinder_viz";
	numFigs = mygetfield( prm, "numFigs", 0 );
	%
	gVals = abs(fVals);
	[ foo, nOfPtWiseMin ] = min( gVals )
	%
	if (isempty(bigSLo))
		bigSLo = xVals(nOfPtWiseMin-1);
	end
	if (isempty(bigSHi))
		bigSHi = xVals(nOfPtWiseMin+1);
	end
	%
	numColors = 1000;
	sizeBigS = 51;
	sizeBigP = 53;
	bigSVals = linspace( bigSLo, bigSHi, sizeBigS );
	bigPVals = linspace( bigPLo, bigPHi, sizeBigP );
	[ bigSMesh, bigPMesh ] = meshgrid( bigSVals, bigPVals );
	tic();
	[ omegaMesh, bigAMesh, bigBMesh, bigCMesh ] = extFinder_getOmegaMesh( bigSMesh, bigPMesh, xVals, fVals, nOfPtWiseMin );
	min(min(omegaMesh))
	max(max(omegaMesh))
	toc();
	%
	numFigs++; figure(numFigs);
	%plot( bigSVals((sizeBigP+1)/2,:), omegaMesh((sizeBigP+1)/2,:) );
	plot( bigSMesh((sizeBigP+[-3,1,3])/2,:)', omegaMesh((sizeBigP+[-3,1,3])/2,:)', 'o-' );
	grid on;
	xlabel( "bigS" );
	ylabel( "omega" );
	title( "omega vs bigS" );
	%
	numFigs++; figure(numFigs);
	plot( bigPMesh(:,(sizeBigS+[-3,1,3])/2), omegaMesh(:,(sizeBigP+[-3,1,3])/2), 'o-' );
	grid on;
	xlabel( "bigP" );
	ylabel( "omega" );
	title( "omega vs bigP" );
	%
	numFigs++; figure(numFigs);
	%contourf( bigSVals, bigPVals, asinh(omegaMesh*1e6)/1e6, 50 );
	%contourf( bigSVals, bigPVals, log(sqrt(eps)*max(max(omegaMesh))+omegaMesh), 50 );
	contourf( bigSVals, bigPVals, (omegaMesh).^0.5, 20 );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "bigS" );
	ylabel( "bigP" );
	title( "omega vs bigS, bigP" );
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( bigSMesh((sizeBigP+[-3,1,3])/2,:)', bigAMesh((sizeBigP+[-3,1,3])/2,:)', 'o-' );
	grid on;
	xlabel( "bigS" );
	ylabel( "bigA" );
	title( "bigA vs bigS" );
	%
	numFigs++; figure(numFigs);
	plot( bigPMesh(:,(sizeBigS+[-3,1,3])/2), bigAMesh(:,(sizeBigP+[-3,1,3])/2), 'o-' );
	grid on;
	xlabel( "bigP" );
	ylabel( "bigA" );
	title( "bigA vs bigP" );
	%
	numFigs++; figure(numFigs);
	contourf( bigSVals, bigPVals, bigAMesh, 20 );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "bigS" );
	ylabel( "bigP" );
	title( "bigA vs bigS, bigP" );
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( bigSMesh((sizeBigP+[-3,1,3])/2,:)', bigBMesh((sizeBigP+[-3,1,3])/2,:)', 'o-' );
	grid on;
	xlabel( "bigS" );
	ylabel( "bigB" );
	title( "bigB vs bigS" );
	%
	numFigs++; figure(numFigs);
	plot( bigPMesh(:,(sizeBigS+[-3,1,3])/2), bigBMesh(:,(sizeBigP+[-3,1,3])/2), 'o-' );
	grid on;
	xlabel( "bigP" );
	ylabel( "bigB" );
	title( "bigB vs bigP" );
	numFigs++; figure(numFigs);
	contourf( bigSVals, bigPVals, bigBMesh, 20 );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "bigS" );
	ylabel( "bigP" );
	title( "bigB vs bigS, bigP" );
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( bigSMesh((sizeBigP+[-3,1,3])/2,:)', bigCMesh((sizeBigP+[-3,1,3])/2,:)', 'o-' );
	grid on;
	xlabel( "bigS" );
	ylabel( "bigC" );
	title( "bigC vs bigS" );
	%
	numFigs++; figure(numFigs);
	plot( bigPMesh(:,(sizeBigS+[-3,1,3])/2), bigCMesh(:,(sizeBigP+[-3,1,3])/2), 'o-' );
	grid on;
	xlabel( "bigP" );
	ylabel( "bigC" );
	title( "bigC vs bigP" );
	numFigs++; figure(numFigs);
	contourf( bigSVals, bigPVals, bigCMesh, 20 );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "bigS" );
	ylabel( "bigP" );
	title( "bigC vs bigS, bigP" );
	%
return;
end
