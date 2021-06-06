%function datOut = extFit_viz( bigX, bigP, rvecX, rvecF, rvecW=[], prm=[] )
	clear;
	setprngstates(26846592); % Massive slowdown.
	%setprngstates(98584480);
	%setprngstates(40121600); % Here, H2 might be better.
	%numPts = 5 + round(abs(randn()*exp(abs(3.0*randn()))))
	numPts = 5 + round(abs(randn()*exp(abs(randn()))))
	bigX_secret = randn()*exp(abs(3.0*randn()))
	bigP_secret = 1.0 + 3.0*abs(randn())
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
	epsX_default = ...
	  min([ sqrt(eps)*(max(rvecX)-min(rvecX)), ...
	   sqrt(sqrt(eps))*min(abs(diff(rvecX))) ]);
	epsP_default = sqrt(eps);
	epsX = mygetfield( prm, "epsX", epsX_default );
	epsP = mygetfield( prm, "epsP", epsP_default );
	%
	omegaDatAtM0 = extFit_calcOmega( bigX0-epsX, bigP0,      rvecX, rvecF, rvecW );
	omegaDatAt0M = extFit_calcOmega( bigX0,      bigP0-epsP, rvecX, rvecF, rvecW );
	omegaDatAt00 = extFit_calcOmega( bigX0,      bigP0,      rvecX, rvecF, rvecW );
	omegaDatAt0P = extFit_calcOmega( bigX0,      bigP0+epsP, rvecX, rvecF, rvecW );
	omegaDatAtP0 = extFit_calcOmega( bigX0+epsX, bigP0,      rvecX, rvecF, rvecW );
	%
	omega0 = omegaDatAt00.omega;
	rvecRhoAtM0 = omegaDatAtM0.rvecRho;
	rvecRhoAt0M = omegaDatAt0M.rvecRho;
	rvecRhoAt00 = omegaDatAt00.rvecRho;
	rvecRhoAt0P = omegaDatAt0P.rvecRho;
	rvecRhoAtP0 = omegaDatAtP0.rvecRho;
	%
	rvecRhoD0 = rvecRhoAt00;
	rvecRhoDX = ( rvecRhoAtP0 - rvecRhoAtM0 ) / ( 2.0*epsX );
	rvecRhoDP = ( rvecRhoAt0P - rvecRhoAt0M ) / ( 2.0*epsP );
	%
	sigma0X = sum( rvecW .* rvecRhoD0 .* rvecRhoDX );
	sigma0P = sum( rvecW .* rvecRhoD0 .* rvecRhoDP );
	sigmaXX = sum( rvecW .* rvecRhoDX .* rvecRhoDX );
	sigmaPP = sum( rvecW .* rvecRhoDP .* rvecRhoDP );
	sigmaXP = sum( rvecW .* rvecRhoDX .* rvecRhoDP );
	%
	vecG = [ sigma0X; sigma0P ];
	matH = [ sigmaXX, sigmaXP; sigmaXP, sigmaPP ];
	matD = diag(abs(diag(matH)));
	matI = eye(2,2);
	%
	%
	%
	sizeBigX = 201;
	sizeBigP = 203;
	rvecBigX = sort( bigX0 + 2.0*(bigX_secret-bigX0)*linspace( -0.5, 1.5, sizeBigX ) );
	rvecBigP = sort( bigP0 + 2.0*(bigP_secret-bigP0)*linspace( -0.5, 1.5, sizeBigP ) );
	[ matBigX, matBigP ] = meshgrid( rvecBigX, rvecBigP );
	%
	dat_calcOmega = extFit_calcOmega_mat( matBigX, matBigP, rvecX, rvecF, rvecW );
	%
	rvecDeltaX = rvecBigX - bigX0;
	rvecDeltaP = rvecBigP - bigP0;
	[ matDeltaX, matDeltaP ] = meshgrid( rvecDeltaX, rvecDeltaP );
	%
	matOmegaModel = omega0 + matH(2,1) * matDeltaX .* matDeltaP ....
	 + matDeltaX .* ( vecG(1) + matH(1,1)*matDeltaX ) ...
	 + matDeltaP .* ( vecG(2) + matH(2,2)*matDeltaP );
	%
	numCurvePts = 53;
	rvecLambda = linspace( 0.0, 1.0, numCurvePts ).^2;
	for n=1:numCurvePts
		lambda = rvecLambda(n);
		if ( lambda < sqrt(eps) )
			matDelta_levenberg(:,n) = -lambda*( matI \ vecG );
		else
			mu = (1.0./lambda) - 1.0;
			matDelta_levenberg(:,n) = -(matH+mu*matI)\vecG;
		end
	end
	%
	%
	numColors = mygetfield( prm, "numColors", 1000 );
	matOmega = dat_calcOmega.matOmega;
	omegaMax = max(max(matOmega));
	omegaMin = 0.0;
	%
	numFigs++; figure(numFigs);
	%contourf( rvecBigX, rvecBigP, matOmega );
	image( rvecBigX, rvecBigP, numColors*(matOmega-omegaMin)/(omegaMax-omegaMin) );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	plot( bigX_secret, bigP_secret, 'w*', 'markersize', 25, 'linewidth', 3 );
	plot( bigX0, bigP0, 'ws', 'markersize', 20, 'linewidth', 2 );
	plot( matDelta_levenberg(1,:)+bigX0, matDelta_levenberg(2,:)+bigP0, 'wo-' );
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omega vs bigX, bigP" );
	%
	numFigs++; figure(numFigs);
	%contourf( rvecBigX, rvecBigP, matOmegaModel1 );
	image( rvecBigX, rvecBigP, numColors*(matOmegaModel-omegaMin)/(omegaMax-omegaMin) );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	plot( bigX_secret, bigP_secret, 'w*', 'markersize', 25, 'linewidth', 3 );
	plot( bigX0, bigP0, 'ws', 'markersize', 20, 'linewidth', 2 );
	plot( matDelta_levenberg(1,:)+bigX0, matDelta_levenberg(2,:)+bigP0, 'wo-' );
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omegaModel vs bigX, bigP" );
return;
%end
