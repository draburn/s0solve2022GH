function extFinder_viz( xVals, fVals, bigSLo, bigSHi, bigPLo, bigPHi, prm=[] )
	commondefs;
	thisFile = "extFinder_viz";
	numFigs = mygetfield( prm, "numFigs", 0 );
	%
	gVals = abs(fVals);
	[ foo, nOfPtWiseMin ] = min( gVals )
	%
	numColors = 1000;
	sizeBigS = 51;
	sizeBigP = 51;
	bigSVals = linspace( bigSLo, bigSHi, sizeBigS );
	bigPVals = linspace( bigPLo, bigPHi, sizeBigP );
	[ bigSMesh, bigPMesh ] = meshgrid( bigSVals, bigPVals );
	tic();
	[ omegaMesh ] = extFinder_getOmegaMesh( bigSMesh, bigPMesh, xVals, fVals, nOfPtWiseMin );
	toc();
	%
	numFigs++; figure(numFigs);
	contourf( bigSVals, bigPVals, ((omegaMesh)), 50 );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	grid on;
	xlabel( "bigS" );
	ylabel( "bigP" );
	title( "omega vs bigS, bigP" );
	%
return;
end
