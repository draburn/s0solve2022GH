	clear;
	commondefs;
	thisFile = "test_extFit_calcOmega";
	numFigs = 0;
	%
	%
	%
	%for n=1:1000 %%%
	for n=1:1 %%%
	setprngstates();
	%setprngstates(6373056); % Very streched.
	%setprngstates(75511728); % A very streched case.
	%setprngstates(13497824); % Blur may help.
	%setprngstates(52235552);
	%setprngstates(91546736); % Blur helps!
	%setprngstates(43730992); % Blur (even 0.1) makes things worse, but basin is broad anyway.
	%setprngstates(4117536); % Blur 0.1 makes things worse.
	%setprngstates(56938368); % Narrow but blur hurts.
	%
	% 1900
	%setprngstates(44883728);  % Rare case where BOTH non-blur and blur are bad.
	%setprngstates(90472720); % Blur is very far.
	%setprngstates(94811216); % Both are bad.
	%setprngstates(28378336); % A below-tol point.
	setprngstates(13835520); % Both are quite bad.
	%
	%
	%
	bigA_secret = randn();
	bigB_secret = randn();
	bigX_secret = randn();
	bigP_secret = 2.0 + abs(randn());
	funchF = @(x)( bigA_secret + bigB_secret * abs( x - bigX_secret ) .^ bigP_secret );
	%
	numPts = round(5 + abs(randn()*exp(randn())));
	xScale = 0.1;
	xVals = sort([ ...
	  bigX_secret-xScale*abs(randn(1,2)), ...
	  bigX_secret+xScale*abs(randn(1,2)), ...
	  bigX_secret+xScale*randn(1,numPts-4) ]);
	fVals = funchF(xVals);
	%
	omegaTol = eps * sum( fVals.^2 );
	%
	%
	%
	sizeBigX = 101;
	sizeBigP = 102;
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
	[ omegaMin, bigPIndexOfMin, bigXIndexOfMin ] = minmin(omegaMesh);
	bigXOfMin = bigXVals(bigXIndexOfMin);
	bigPOfMin = bigPVals(bigPIndexOfMin);
	singlePt_omega = extFit_calcOmega( xVals, fVals, bigXOfMin, bigPOfMin );
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
	bigXOfBlurMin = bigXVals(bigXIndexOfBlurMin+1);
	bigPOfBlurMin = bigPVals(bigPIndexOfBlurMin+1);
	%
	%
	bigXRes = abs(bigX_secret-bigXOfMin);
	bigPRes = abs(bigP_secret-bigPOfMin);
	bigXBlurRes = abs(bigX_secret-bigXOfBlurMin);
	bigPBlurRes = abs(bigP_secret-bigPOfBlurMin);
	if (  (bigXRes > 0.2 || bigPRes > 1.0)  &&  (bigXBlurRes > 0.2 || bigPBlurRes > 1.0)  )
		msg( thisFile, __LINE__, "BREAK!" );
		break;
	else
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
