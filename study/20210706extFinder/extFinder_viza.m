function extFinder_viza( xVals, fVals, bigS, bigP, prm=[] )
	commondefs;
	thisFile = "extFinder_viza";
	numFigs = mygetfield( prm, "numFigs", 0 );
	%
	gVals = abs(fVals);
	[ foo, nOfPtWiseMin ] = min( gVals )
	n = nOfPtWiseMin;
	bigDelta = xVals(n+1) - xVals(n-1);
	bigG0 = gVals(n);
	bigG1 = gVals(n+1) + gVals(n-1) - 2.0*gVals(n);
	%
	[ bigA, bigB, bigC ] = extFinder_getFit( bigS, bigP, xVals, gVals, nOfPtWiseMin, prm )
	funchGModel = @(x)( bigG0 + bigG1*( ...
	  bigA*abs((x-bigS)/bigDelta).^bigP + bigB*(x-bigS)/bigDelta + bigC ) );
	semilogy( ...
	  xVals, gVals, 'o-', ...
	  xVals, funchGModel(xVals), 'x-' );
	grid on;
	xlabel( "x" );
	ylabel( "|F|" );
	title( "|F| vs x" );
	%
return;
end
