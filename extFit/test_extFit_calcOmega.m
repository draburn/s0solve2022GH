	clear;
	commondefs;
	thisFile = "test_extFit_calcOmega";
	numFigs = 0;
	%
	%
	%
	%for n=1:1000 %%%
	%for n=1:1 %%%
	setprngstates();
	%setprngstates(6373056); % Very streched.
	%setprngstates(75511728); % A very streched case.
	%setprngstates(13497824); % Blur may help.
	%setprngstates(52235552);
	%setprngstates(91546736); % Blur helps!
	%setprngstates(43730992); % Blur (even 0.1) makes things worse, but basin is broad anyway.
	%setprngstates(4117536); % Blur 0.1 makes things worse.
	setprngstates(56938368); % Narrow but blur hurts.
	%
	% 1900
	%setprngstates(44883728);  % Rare case where BOTH non-blur and blur are bad.
	%setprngstates(90472720); % Blur is very far.
	%setprngstates(94811216); % Both are bad.
	%setprngstates(28378336); % A below-tol point.
	%setprngstates(13835520); % Both are quite bad.
	%
	%
	%
	bigA_secret = randn()
	bigB_secret = randn()
	bigX_secret = randn()
	bigP_secret = 2.0 + abs(randn())
	funchF = @(x)( bigA_secret + bigB_secret * abs( x - bigX_secret ) .^ bigP_secret );
	%
	numPts = round(5 + abs(randn()*exp(randn())));
	xVals = sort([ ...
	  bigX_secret-abs(randn(1,2)), ...
	  bigX_secret+abs(randn(1,2)), ...
	  bigX_secret+randn(1,numPts-4) ]);
	fVals = funchF(xVals);
	%
	wVals = ones(size(xVals));
	%wVals = 1.0./sqrt(sqrt( eps*max(abs(fVals)) + abs(fVals) ));
	%
	omegaTol = eps * sum( fVals.^2 );
	%
	%
	%
	sizeBigX = 501;
	sizeBigP = 502;
	indexOfExt = 2;
	while ( (fVals(indexOfExt+1)-fVals(indexOfExt)) * (fVals(2)-fVals(1)) > 0 )
		indexOfExt++;
	end
	bigXLo = xVals(indexOfExt-1);
	bigXHi = xVals(indexOfExt+1);
	%bigXLo = min(xVals);
	%bigXHi = max(xVals);
	bigPLo = 0.5;
	bigPHi = 10.0;
	bigXVals = linspace( bigXLo, bigXHi, sizeBigX );
	bigPVals = linspace( bigPLo, bigPHi, sizeBigP );
	[ bigXMesh, bigPMesh ] = meshgrid( bigXVals, bigPVals );
	%
	numColors = 1000;
	%
	%
	%
	[ omegaMesh, rhoAry3, bigAMesh, bigBMesh ] = extFit_calcOmega( ...
	  xVals, fVals, bigXMesh, bigPMesh, wVals );
	%
	belowTolMesh = (omegaMesh<omegaTol);
	if ( sum(sum(belowTolMesh)) > 0 )
		% Pretty unlikely.
		% But, in this case, don't blur.
		% And, don't bother to check for multiples.
		msg( thisFile, __LINE__, "Found a below-tol point!" );
	end
	%
	[ omegaMin, bigPIndexOfMin, bigXIndexOfMin ] = minmin(omegaMesh);
	bigXOfMin = bigXVals(bigXIndexOfMin)
	bigPOfMin = bigPVals(bigPIndexOfMin)
	[ ofMin_omega, ofMin_rhoVec, ofMin_bigA, ofMin_bigB ] = extFit_calcOmega( ...
	  xVals, fVals, bigXOfMin, bigPOfMin, wVals );
	echo__ofMin_omega = ofMin_omega
	echo__ofMin_bigA = ofMin_bigA
	echo__ofMin_bigB = ofMin_bigB
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
	bigXRes = abs(bigX_secret-bigXOfMin);
	bigPRes = abs(bigP_secret-bigPOfMin);
	bigXBlurRes = abs(bigX_secret-bigXOfBlurMin);
	bigPBlurRes = abs(bigP_secret-bigPOfBlurMin);
	%if (  (bigXRes > 0.2 || bigPRes > 1.0)  &&  (bigXBlurRes > 0.2 || bigPBlurRes > 1.0)  )
	%	msg( thisFile, __LINE__, "BREAK!" );
	%	break;
	%else
	%	continue;
	%end
	%end %%%
	%
	numFigs++; figure(numFigs);
	%imagesc( bigXVals, bigPVals, omegaMesh );
	imagesc( bigXVals(2:end-1), bigPVals(2:end-1), omegaBlurMesh );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	plot( bigX_secret, bigP_secret, 'wo', 'linewidth', 2, 'markersize', 25 );
	plot( bigXOfMin, bigPOfMin, 'wx', 'linewidth', 2, 'markersize', 25 );
	plot( bigXOfBlurMin, bigPOfBlurMin, 'w+', 'linewidth', 2, 'markersize', 25 );
	doActualSolve = true;
	if (doActualSolve)
		%bigX0 = ( xVals(indexOfExt-1) + xVals(indexOfExt+1) )/2.0;
		%bigP0 = 2.0;
		bigX0 = bigXOfBlurMin;
		bigP0 = bigPOfBlurMin;
		prm_extFit.omegaTol = omegaTol;
		[ dat_extFit, retCode ] = extFit( bigX0, bigP0, xVals, fVals, [], prm_extFit, wVals );
		plot( dat_extFit.bigX, dat_extFit.bigP, 'r^', 'linewidth', 2, 'markersize', 15 );
	end
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omega vs bigX, bigP" );
	%
	numFigs++; figure(numFigs);
	viz_xVals = linspace( min(xVals), max(xVals), 1000 );
	viz_fVals = funchF( viz_xVals );
	plot( ...
	  xVals, fVals, 'kx', 'linewidth', 4, 'markersize', 15, ...
	  viz_xVals, viz_fVals, 'o-', 'linewidth', 1, 'markersize', 4 );
	hold on;
	if (doActualSolve)
		solver_bigX = dat_extFit.bigX;
		solver_bigP = dat_extFit.bigP;
		[ solver_omega, solver_rho, solver_bigA, solver_bigB ] = extFit_calcOmega( ...
		  xVals, fVals, solver_bigX, solver_bigP, wVals );
		funchFModel = @(x)( solver_bigA + solver_bigB * abs( x - solver_bigX ).^solver_bigP );
		viz_fModelVals = funchFModel( viz_xVals );
		plot( viz_xVals, viz_fModelVals, 'g^-', 'linewidth', 1, 'markersize', 4 );
	end
	hold off;
	grid on;
	xlabel( 'x' );
	ylabel( 'f' );
	title( 'points to fit' );
