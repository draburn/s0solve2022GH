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
