function extFit_viz( xVals, fVals, bigXLo, bigXHi, bigPLo, bigPHi, prm=[] )
	commondefs;
	thisFile = "extFit_viz";
	numFigs = mygetfield( prm, "numFigs", 0 );
	%
	omegaTol = eps * sum( fVals.^2 );
	%
	numColors = 1000;
	sizeBigX = 51;
	sizeBigP = 52;
	bigXVals = linspace( bigXLo, bigXHi, sizeBigX );
	bigPVals = linspace( bigPLo, bigPHi, sizeBigP );
	[ bigXMesh, bigPMesh ] = meshgrid( bigXVals, bigPVals );
	wVals = mygetfield( prm, "wVals", [] );
	%msg( thisFile, __LINE__, "" );
	%msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
	%extFit_calcOmega( xVals, fVals, 1.0*ones(2,2), 2.1*ones(2,2), wVals )
	[ omegaMesh ] = extFit_calcOmega( xVals, fVals, bigXMesh, bigPMesh, wVals );
	%msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
	%msg( thisFile, __LINE__, "" );
	%
	belowTolMesh = (omegaMesh<omegaTol);
	if ( sum(sum(belowTolMesh)) > 0 )
		% Pretty unlikely.
		% But, in this case, don't blur.
		% And, don't bother to check for multiples.
		msg( thisFile, __LINE__, "Found a below-tol point!" );
	end
	%
	omegaMesh = log(eps+omegaMesh);
	%
	[ omegaMin, bigPIndexOfMin, bigXIndexOfMin ] = minmin(omegaMesh);
	bigXOfMin = bigXVals(bigXIndexOfMin)
	bigPOfMin = bigPVals(bigPIndexOfMin)
	[ ofMin_omega, ofMin_rhoVec, ofMin_bigA, ofMin_bigB ] = extFit_calcOmega( xVals, fVals, bigXOfMin, bigPOfMin );
	%
	blurCoeff = 0.5;
	tempMesh = ...
	   (1.0-blurCoeff) * omegaMesh(2:end-1,:)  ...
	 + (blurCoeff/2.0) * ( omegaMesh(3:end,:) + omegaMesh(1:end-2,:) );
	omegaBlurMesh = ...
	   (1.0-blurCoeff) * tempMesh(:,2:end-1)  ...
	 + (blurCoeff/2.0) * ( tempMesh(:,3:end) + tempMesh(:,1:end-2) );
	%
	[ omegaBlurMin, bigPIndexOfBlurMin, bigXIndexOfBlurMin ] = minmin(omegaBlurMesh);
	bigXOfBlurMin = bigXVals(bigXIndexOfBlurMin+1)
	bigPOfBlurMin = bigPVals(bigPIndexOfBlurMin+1)
	%
	%
	%bigXRes = abs(bigX_secret-bigXOfMin);
	%bigPRes = abs(bigP_secret-bigPOfMin);
	%bigXBlurRes = abs(bigX_secret-bigXOfBlurMin);
	%bigPBlurRes = abs(bigP_secret-bigPOfBlurMin);
	%%%if (  (bigXRes > 0.2 || bigPRes > 1.0)  &&  (bigXBlurRes > 0.2 || bigPBlurRes > 1.0)  )
	%%%	msg( thisFile, __LINE__, "BREAK!" );
	%%%	break;
	%%%else
	%%%	continue;
	%%%end
	%
	numFigs++; figure(numFigs);
	imagesc( bigXVals, bigPVals, omegaMesh );
	%imagesc( bigXVals(2:end-1), bigPVals(2:end-1), omegaBlurMesh );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	%plot( bigX_secret, bigP_secret, 'wo', 'linewidth', 2, 'markersize', 25 );
	plot( bigXOfMin, bigPOfMin, 'wx', 'linewidth', 2, 'markersize', 25 );
	plot( bigXOfBlurMin, bigPOfBlurMin, 'w+', 'linewidth', 2, 'markersize', 25 );
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omega vs bigX, bigP" );
	%
return;
end
