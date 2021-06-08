	clear;
	commondefs;
	thisFile = "test_extFit_calcOmega";
	numFigs = 0;
	%
	%
	%
	for n=1:1000 %%%
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
	omegaTol = eps * sum( fVals.^2 );
	%
	%
	%
	sizeBigX = 101;
	sizeBigP = 102;
	bigXLo = min(xVals);
	bigXHi = max(xVals);
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
	[ omegaMesh, rhoAry3, bigAMesh, bigBMesh ] = extFit_calcOmega( xVals, fVals, bigXMesh, bigPMesh );
	%
	belowTolMesh = (omegaMesh<omegaTol);
	if ( sum(sum(belowTolMesh)) > 0 )
		% Pretty unlikely.
		% But, in this case, don't blur.
		% And, don't bother to check for multiples.
		msg( thisFile, __LINE__, "Found a below-tol point!" );
	end
	%
	[ omegaMin, bigPIndexOfMin, bigXIndexOfMin ] = minmin(omegaMesh)
	bigXOfMin = bigXVals(bigXIndexOfMin)
	bigPOfMin = bigPVals(bigPIndexOfMin)
	singlePt_omega = extFit_calcOmega( xVals, fVals, bigXOfMin, bigPOfMin )
	%
	blurCoeff = 0.01;
	tempMesh = (  (1.0-blurCoeff)*omegaMesh(2:end-1,:)  +  blurCoeff * ...
	  ( omegaMesh(3:end,:) + omegaMesh(1:end-2,:) )  ) / 3.0;
	omegaBlurMesh = (  (1.0-blurCoeff)*tempMesh(:,2:end-1)  +  blurCoeff * ...
	  ( tempMesh(:,3:end)  + tempMesh(:,1:end-2)  )  ) / 3.0;
	%
	[ omegaBlurMin, bigPIndexOfBlurMin, bigXIndexOfBlurMin ] = minmin(omegaBlurMesh)
	bigXOfBlurMin = bigXVals(bigXIndexOfBlurMin+1)
	bigPOfBlurMin = bigPVals(bigPIndexOfBlurMin+1)
	%
	%
	bigXRes = abs(bigX_secret-bigXOfBlurMin);
	bigPRes = abs(bigP_secret-bigPOfBlurMin);
	if ( bigXRes > 0.2 || bigPRes > 1.0 )
		"BREAK!"
		break;
	else
		"CONTINUE!"
		continue;
	end
	end %%%
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
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omega vs bigX, bigP" );
