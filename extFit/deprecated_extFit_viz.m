%function datOut = extFit_viz( bigX, bigP, rvecX, rvecF, rvecW=[], prm=[] )
	clear;
	setprngstates(98584480);
	%setprngstates(40121600); % Here, H2 might be better.
	numPts = 5 + round(abs(randn()*exp(abs(3.0*randn()))))
	bigX_secret = randn()*exp(abs(3.0*randn()))
	bigP_secret = 1.0 + 3.0*abs(randn());
	bigA_secret = randn()*exp(abs(3.0*randn()));
	bigB_secret = randn()*exp(abs(3.0*randn()));
	rvecX = bigX_secret + randn(1,numPts);
	funchF = @(x)( bigA_secret + bigB_secret * abs( rvecX - bigX_secret ).^bigP_secret );
	rvecF = funchF(rvecX);
	rvecW = [];
	prm = [];
	bigX0 = bigX_secret + randn()
	bigP0 = 2.0
	%
	commondefs;
	thisFile = "extFit_viz";
	numFigs = 0;
	%
	msg( thisFile, __LINE__, "This makes use of H2, largerly illustrating why I should abandon it." );
	%
	if (isempty(rvecW))
		rvecW = ones(size(rvecX));
	end
	numPts = size(rvecX,2);
	assert( isrealscalar(bigX0) );
	assert( isrealscalar(bigP0) );
	assert( isrealarray(rvecX,[1,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(rvecW,[1,numPts]) );
	%
	%
	%
	prm_calcGradHess = mygetfield( prm, "prm_calcGradHess", [] );
	dat_calcGradHess = extFit_calcGradHess( bigX0, bigP0, rvecX, rvecF, rvecW, prm_calcGradHess );
	omega0 = dat_calcGradHess.omega0;
	vecG = dat_calcGradHess.vecG;
	matH1 = dat_calcGradHess.matH1;
	matH2 = dat_calcGradHess.matH2;
	matD1 = diag(abs(diag(matH1)));
	matD2 = diag(abs(diag(matH2)));
	matI = eye(2,2);
	%
	%
	%
	sizeBigX = 201;
	sizeBigP = 203;
	deltaX_secret = abs( bigX0 - bigX_secret );
	deltaP_secret = abs( bigP0 - bigP_secret );
	rvecBigX = bigX0 + 3.0*deltaX_secret*linspace( -1.0, 1.0, sizeBigX );
	rvecBigP = bigP0 + 30.0*deltaP_secret*linspace( -1.0, 1.0, sizeBigP );
	[ matBigX, matBigP ] = meshgrid( rvecBigX, rvecBigP );
	%
	dat_calcOmega = extFit_calcOmega_mat( matBigX, matBigP, rvecX, rvecF, rvecW );
	%
	rvecDeltaX = rvecBigX - bigX0;
	rvecDeltaP = rvecBigP - bigP0;
	[ matDeltaX, matDeltaP ] = meshgrid( rvecDeltaX, rvecDeltaP );
	%
	matOmegaModel1 = omega0 + matH1(2,1) * matDeltaX .* matDeltaP ....
	 + matDeltaX .* ( vecG(1) + matH1(1,1)*matDeltaX ) ...
	 + matDeltaP .* ( vecG(2) + matH1(2,2)*matDeltaP );
	matOmegaModel2 = omega0 + matH2(2,1) * matDeltaX .* matDeltaP ....
	 + matDeltaX .* ( vecG(1) + matH2(1,1)*matDeltaX ) ...
	 + matDeltaP .* ( vecG(2) + matH2(2,2)*matDeltaP );
	%
	prm_findMuOfOmega = mygetfield( prm, "prm_findMuOfOmega", [] );
	muMin_h2_d1 = extFit_findMuOfOmega( 0.0, omega0, vecG, matH2, matD1, prm_findMuOfOmega )
	muMin_h2_d2 = extFit_findMuOfOmega( 0.0, omega0, vecG, matH2, matD2, prm_findMuOfOmega )
	muMin_h2_i  = extFit_findMuOfOmega( 0.0, omega0, vecG, matH2, matI, prm_findMuOfOmega )
	%
	numCurvePts = 53;
	rvecLambda = linspace( 0.0, 1.0, numCurvePts ).^2;
	for n=1:numCurvePts
		lambda = rvecLambda(n);
		if ( lambda < sqrt(eps) )
			matDelta_h1_i(:,n) = -lambda*( matI \ vecG );
			matDelta_h2_i(:,n) = -lambda*( matI \ vecG );
		else
			mu = (1.0./lambda) - 1.0;
			matDelta_h1_i(:,n) = -(matH1+mu*matI)\vecG;
			matDelta_h2_i(:,n) = -(matH2+(mu+muMin_h2_i)*matI)\vecG;
			matDelta_h2_d1(:,n) = -(matH2+(mu+muMin_h2_d1)*matD1)\vecG;
			matDelta_h2_d2(:,n) = -(matH2+(mu+muMin_h2_d2)*matD2)\vecG;
		end
	end
	%
	%
	numColors = mygetfield( prm, "numColors", 1000 );
	matOmega = dat_calcOmega.matOmega;
	omegaMax = max(max(matOmega));
	omegaMin = 0.0;
	%omegaMax = max(max(matOmegaModel2));
	%omegaMin = min(min(matOmegaModel2));
	%
	numFigs++; figure(numFigs);
	%contourf( rvecBigX, rvecBigP, matOmega );
	image( rvecBigX, rvecBigP, numColors*(matOmega-omegaMin)/(omegaMax-omegaMin) );
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	plot( bigX_secret, bigP_secret, 'w*', 'markersize', 25, 'linewidth', 3 );
	plot( bigX0, bigP0, 'ws', 'markersize', 20, 'linewidth', 2 );
	plot( matDelta_h1_i(1,:)+bigX0, matDelta_h1_i(2,:)+bigP0, 'wo-' );
	plot( matDelta_h2_i(1,:)+bigX0, matDelta_h2_i(2,:)+bigP0, 'ko-' );
	plot( matDelta_h2_d1(1,:)+bigX0, matDelta_h2_d1(2,:)+bigP0, 'ro-' );
	plot( matDelta_h2_d2(1,:)+bigX0, matDelta_h2_d2(2,:)+bigP0, 'bo-' );
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omega vs bigX, bigP" );
	%
	numFigs++; figure(numFigs);
	%contourf( rvecBigX, rvecBigP, matOmegaModel1 );
	image( rvecBigX, rvecBigP, numColors*(matOmegaModel1-omegaMin)/(omegaMax-omegaMin) );
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	plot( bigX_secret, bigP_secret, 'w*', 'markersize', 25, 'linewidth', 3 );
	plot( bigX0, bigP0, 'ws', 'markersize', 20, 'linewidth', 2 );
	plot( matDelta_h1_i(1,:)+bigX0, matDelta_h1_i(2,:)+bigP0, 'wo-' );
	plot( matDelta_h2_i(1,:)+bigX0, matDelta_h2_i(2,:)+bigP0, 'ko-' );
	plot( matDelta_h2_d1(1,:)+bigX0, matDelta_h2_d1(2,:)+bigP0, 'ro-' );
	plot( matDelta_h2_d2(1,:)+bigX0, matDelta_h2_d2(2,:)+bigP0, 'bo-' );
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omegaModel1 vs bigX, bigP" );
	%
	numFigs++; figure(numFigs);
	%contourf( rvecBigX, rvecBigP, matOmegaModel2 );
	image( rvecBigX, rvecBigP, numColors*(matOmegaModel2-omegaMin)/(omegaMax-omegaMin) );
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	plot( bigX_secret, bigP_secret, 'w*', 'markersize', 25, 'linewidth', 3 );
	plot( bigX0, bigP0, 'ws', 'markersize', 20, 'linewidth', 2 );
	plot( matDelta_h1_i(1,:)+bigX0, matDelta_h1_i(2,:)+bigP0, 'wo-' );
	plot( matDelta_h2_i(1,:)+bigX0, matDelta_h2_i(2,:)+bigP0, 'ko-' );
	plot( matDelta_h2_d1(1,:)+bigX0, matDelta_h2_d1(2,:)+bigP0, 'ro-' );
	plot( matDelta_h2_d2(1,:)+bigX0, matDelta_h2_d2(2,:)+bigP0, 'bo-' );
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omegaModel2 vs bigX, bigP" );
	%
return;
%end
