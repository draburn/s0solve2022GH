%function datOut = extFit_viz( bigX, bigP, rvecX, rvecF, rvecW=[], prm=[] )
	clear;
	setprngstates();
	%setprngstates(10743664); % Image plots glitched. Turn off axis(suqare), "a", zoom out.
	%setprngstates(77077008);
	%setprngstates(5932672)
	%setprngstates(5932672);
	%setprngstates(30660864); % Start on wrong side of a point???
	%setprngstates(58467872);
	%setprngstates(26846592); % Massive slowdown. And, ?!?!?!?!
	%setprngstates(98584480);
	%setprngstates(40121600); % Here, H2 might be better.
	%numPts = 5 + round(abs(randn()*exp(abs(3.0*randn()))))
	numPts = 5 + round(abs(randn()*exp(abs(randn()))))
	bigX_secret = randn()*exp(abs(3.0*randn()))
	bigP_secret = 1.0 + 3.0*abs(randn())
	bigA_secret = randn()*exp(abs(3.0*randn()));
	bigB_secret = randn()*exp(abs(3.0*randn()));
	rvecX = sort(bigX_secret + randn(1,numPts));
	funchF = @(x)( bigA_secret + bigB_secret * abs( x - bigX_secret ).^bigP_secret );
	rvecF = funchF(rvecX);
	rvecW = [];
	prm = [];
	index0 = 1;
	if ( bigB_secret > 0 )
	while (1)
		if ( (index0==numPts) )
			break;
		elseif ( rvecF(index0+1) > rvecF(index0) )
			break;
		else
			index0++;
			continue;
		end
	end
	elseif ( bigB_secret < 0 )
	while (1)
		if ( (index0==numPts) )
			break;
		elseif ( rvecF(index0+1) < rvecF(index0) )
			break;
		else
			index0++;
			continue;
		end
	end
	end
	if ( 1==index0 )
		bigX0 = (rvecX(1)+rvecX(2))/2.0
	elseif ( numPts==index0)
		bigX0 = (rvecX(numPts)+rvecX(numPts-1))/2.0
	else
		bigX0 = (rvecX(index0+1)+rvecX(index0-1))/2.0
	end
	bigP0 = 2.0
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
	rvecX_first = linspace(min(rvecX),max(rvecX),101);
	rvecF_first = funchF(rvecX_first);
	numFigs++; figure(numFigs);
	plot( ...
	  rvecX, rvecF, 'o', 'markersize', 20, ...
	  rvecX_first, rvecF_first, '*-', ...
	  bigX0*[1,1], [min(rvecF_first),max(rvecF_first)], 'k-' );
	xlabel( "x" );
	ylabel( "f" );
	title( "f vx x" );
	grid on;
	%
	%
	%
	prm_calcGradHess = mygetfield( prm, "prm_calcGradHess", [] );
	dat_calcGradHess = extFit_calcGradHess( bigX0, bigP0, rvecX, rvecF, rvecW, prm_calcGradHess );
	omega0 = dat_calcGradHess.omega0;
	vecG = dat_calcGradHess.vecG;
	matH1 = dat_calcGradHess.matH1;
	assert( 0==sum( diag(matH1)<=0.0 ) );
	matH2 = dat_calcGradHess.matH2;
	matD0 = eye(2,2);
	matD1 = diag(diag(matH1));
	matD2 = diag(abs(diag(matH2)));
	%
	prm_findMuOfOmega = mygetfield( prm, "prm_findMuOfOmega", [] );
	muMin_h1d0 = 0.0;
	muMin_h1d1 = 0.0;
	[ muMin_h2d0, retCode_d0 ] = extFit_findMuOfOmega( 0.0, omega0, vecG, matH2, matD0, prm_findMuOfOmega );
	[ muMin_h2d1, retCode_d1 ] = extFit_findMuOfOmega( 0.0, omega0, vecG, matH2, matD1, prm_findMuOfOmega );
	[ muMin_h2d2, retCode_d2 ] = extFit_findMuOfOmega( 0.0, omega0, vecG, matH2, matD2, prm_findMuOfOmega );
	% DRaburn 2021.06.06.
	% We could reject results based on retCode.
	% But, we probably don't want to ultimately use H2 anyway.
	% So, this is just for viz.
	%
	muCoeff_h1d0 = norm( matD0 \ vecG ) / norm( (matH1+(muMin_h1d0*matD0))\vecG );
	muCoeff_h1d1 = norm( matD1 \ vecG ) / norm( (matH1+(muMin_h1d1*matD1))\vecG );
	muCoeff_h2d0 = norm( matD0 \ vecG ) / norm( (matH2+(muMin_h2d0*matD0))\vecG );
	muCoeff_h2d1 = norm( matD1 \ vecG ) / norm( (matH2+(muMin_h2d1*matD1))\vecG );
	muCoeff_h2d2 = norm( matD2 \ vecG ) / norm( (matH2+(muMin_h2d2*matD2))\vecG );
	%
	%
	numCurvePts = 101;
	rvecLambda = linspace( 0.0, 1.0, numCurvePts ).^0.5;
	rvecMu_h1d0 = muMin_h1d0 + muCoeff_h1d0 * ( (1.0./rvecLambda) - 1.0 );
	rvecMu_h1d1 = muMin_h1d1 + muCoeff_h1d1 * ( (1.0./rvecLambda) - 1.0 );
	rvecMu_h2d0 = muMin_h2d0 + muCoeff_h2d0 * ( (1.0./rvecLambda) - 1.0 );
	rvecMu_h2d1 = muMin_h2d1 + muCoeff_h2d1 * ( (1.0./rvecLambda) - 1.0 );
	rvecMu_h2d2 = muMin_h2d2 + muCoeff_h2d2 * ( (1.0./rvecLambda) - 1.0 );
	for n=1:numCurvePts
		lambda = rvecLambda(n);
		if ( lambda < sqrt(eps) )
			matDelta_h1d0(:,n) = -lambda * ( matD0 \ vecG ) / muCoeff_h1d0;
			matDelta_h1d1(:,n) = -lambda * ( matD1 \ vecG ) / muCoeff_h1d1;
			matDelta_h2d0(:,n) = -lambda * ( matD0 \ vecG ) / muCoeff_h2d0;
			matDelta_h2d1(:,n) = -lambda * ( matD1 \ vecG ) / muCoeff_h2d1;
			matDelta_h2d2(:,n) = -lambda * ( matD2 \ vecG ) / muCoeff_h2d2;
		else
			matDelta_h1d0(:,n) = -( matH1 + rvecMu_h1d0(n)*matD0 ) \ vecG;
			matDelta_h1d1(:,n) = -( matH1 + rvecMu_h1d1(n)*matD1 ) \ vecG;
			matDelta_h2d0(:,n) = -( matH2 + rvecMu_h2d0(n)*matD0 ) \ vecG;
			matDelta_h2d1(:,n) = -( matH2 + rvecMu_h2d1(n)*matD1 ) \ vecG;
			matDelta_h2d2(:,n) = -( matH2 + rvecMu_h2d2(n)*matD2 ) \ vecG;
		end
		%
		rvecDeltaNorm_h1d0(n) = norm(matDelta_h1d0(:,n));
		rvecDeltaNorm_h1d1(n) = norm(matDelta_h1d1(:,n));
		rvecDeltaNorm_h2d0(n) = norm(matDelta_h2d0(:,n));
		rvecDeltaNorm_h2d1(n) = norm(matDelta_h2d1(:,n));
		rvecDeltaNorm_h2d2(n) = norm(matDelta_h2d2(:,n));
		%
		rvecOmegaM1_h1d0(n) = omega0 + vecG'*matDelta_h1d0(:,n) + 0.5*matDelta_h1d0(:,n)'*matH1*matDelta_h1d0(:,n);
		rvecOmegaM1_h1d1(n) = omega0 + vecG'*matDelta_h1d1(:,n) + 0.5*matDelta_h1d1(:,n)'*matH1*matDelta_h1d1(:,n);
		rvecOmegaM1_h2d0(n) = omega0 + vecG'*matDelta_h2d0(:,n) + 0.5*matDelta_h2d0(:,n)'*matH1*matDelta_h2d0(:,n);
		rvecOmegaM1_h2d1(n) = omega0 + vecG'*matDelta_h2d1(:,n) + 0.5*matDelta_h2d1(:,n)'*matH1*matDelta_h2d1(:,n);
		rvecOmegaM1_h2d2(n) = omega0 + vecG'*matDelta_h2d2(:,n) + 0.5*matDelta_h2d2(:,n)'*matH1*matDelta_h2d2(:,n);
		%
		rvecOmegaM2_h1d0(n) = omega0 + vecG'*matDelta_h1d0(:,n) + 0.5*matDelta_h1d0(:,n)'*matH2*matDelta_h1d0(:,n);
		rvecOmegaM2_h1d1(n) = omega0 + vecG'*matDelta_h1d1(:,n) + 0.5*matDelta_h1d1(:,n)'*matH2*matDelta_h1d1(:,n);
		rvecOmegaM2_h2d0(n) = omega0 + vecG'*matDelta_h2d0(:,n) + 0.5*matDelta_h2d0(:,n)'*matH2*matDelta_h2d0(:,n);
		rvecOmegaM2_h2d1(n) = omega0 + vecG'*matDelta_h2d1(:,n) + 0.5*matDelta_h2d1(:,n)'*matH2*matDelta_h2d1(:,n);
		rvecOmegaM2_h2d2(n) = omega0 + vecG'*matDelta_h2d2(:,n) + 0.5*matDelta_h2d2(:,n)'*matH2*matDelta_h2d2(:,n);
	end
	%
	dat_calcOmega_h1d0 = extFit_calcOmega_mat( ...
	  bigX0 + matDelta_h1d0(1,:), ...
	  bigP0 + matDelta_h1d0(2,:), ...
	  rvecX, rvecF, rvecW );
	dat_calcOmega_h1d1 = extFit_calcOmega_mat( ...
	  bigX0 + matDelta_h1d1(1,:), ...
	  bigP0 + matDelta_h1d1(2,:), ...
	  rvecX, rvecF, rvecW );
	dat_calcOmega_h2d0 = extFit_calcOmega_mat( ...
	  bigX0 + matDelta_h2d0(1,:), ...
	  bigP0 + matDelta_h2d0(2,:), ...
	  rvecX, rvecF, rvecW );
	dat_calcOmega_h2d1 = extFit_calcOmega_mat( ...
	  bigX0 + matDelta_h2d1(1,:), ...
	  bigP0 + matDelta_h2d1(2,:), ...
	  rvecX, rvecF, rvecW );
	dat_calcOmega_h2d2 = extFit_calcOmega_mat( ...
	  bigX0 + matDelta_h2d2(1,:), ...
	  bigP0 + matDelta_h2d2(2,:), ...
	  rvecX, rvecF, rvecW );
	%
	rvecOmegaAc_h1d0 = dat_calcOmega_h1d0.matOmega;
	rvecOmegaAc_h1d1 = dat_calcOmega_h1d1.matOmega;
	rvecOmegaAc_h2d0 = dat_calcOmega_h2d0.matOmega;
	rvecOmegaAc_h2d1 = dat_calcOmega_h2d1.matOmega;
	rvecOmegaAc_h2d2 = dat_calcOmega_h2d2.matOmega;
	%
	%
	%
	plot_color_h1d0 = [ 0.7, 0.0, 0.0 ];
	plot_color_h1d1 = [ 0.6, 0.5, 0.0 ];
	plot_color_h2d0 = [ 0.0, 0.6, 0.0 ];
	plot_color_h2d1 = [ 0.0, 0.0, 0.8 ];
	plot_color_h2d2 = [ 0.6, 0.0, 0.7 ];
	plot_marker_h1d0 = 'o-';
	plot_marker_h1d1 = 's-';
	plot_marker_h2d0 = 'x-';
	plot_marker_h2d1 = '+-';
	plot_marker_h2d2 = '*-';
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rvecLambda, rvecDeltaNorm_h1d0, plot_marker_h1d0, 'color', plot_color_h1d0, ...
	  rvecLambda, rvecDeltaNorm_h1d1, plot_marker_h1d1, 'color', plot_color_h1d1, ...
	  rvecLambda, rvecDeltaNorm_h2d0, plot_marker_h2d0, 'color', plot_color_h2d0, ...
	  rvecLambda, rvecDeltaNorm_h2d1, plot_marker_h2d1, 'color', plot_color_h2d1, ...
	  rvecLambda, rvecDeltaNorm_h2d2, plot_marker_h2d2, 'color', plot_color_h2d2 );
	legend( ...
	  "H1D0", ...
	  "H1D1", ...
	  "H2D0", ...
	  "H2D1", ...
	  "H2D2", ...
	  "location", "north" );
	xlabel( "lambda" );
	ylabel( "delta norm" );
	title( "delta norm vs lambda" );
	grid on;
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rvecDeltaNorm_h1d0, rvecOmegaAc_h1d0, plot_marker_h1d0, 'color', plot_color_h1d0, ...
	  rvecDeltaNorm_h1d1, rvecOmegaAc_h1d1, plot_marker_h1d1, 'color', plot_color_h1d1, ...
	  rvecDeltaNorm_h2d0, rvecOmegaAc_h2d0, plot_marker_h2d0, 'color', plot_color_h2d0, ...
	  rvecDeltaNorm_h2d1, rvecOmegaAc_h2d1, plot_marker_h2d1, 'color', plot_color_h2d1, ...
	  rvecDeltaNorm_h2d2, rvecOmegaAc_h2d2, plot_marker_h2d2, 'color', plot_color_h2d2 );
	legend( ...
	  "H1D0", ...
	  "H1D1", ...
	  "H2D0", ...
	  "H2D1", ...
	  "H2D2", ...
	  "location", "north" );
	xlabel( "delta norm" );
	ylabel( "omega actual" );
	title( "omega actual vs delta norm" );
	grid on;
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rvecDeltaNorm_h1d0, rvecOmegaM1_h1d0, plot_marker_h1d0, 'color', plot_color_h1d0, ...
	  rvecDeltaNorm_h1d1, rvecOmegaM1_h1d1, plot_marker_h1d1, 'color', plot_color_h1d1, ...
	  rvecDeltaNorm_h2d0, rvecOmegaM1_h2d0, plot_marker_h2d0, 'color', plot_color_h2d0, ...
	  rvecDeltaNorm_h2d1, rvecOmegaM1_h2d1, plot_marker_h2d1, 'color', plot_color_h2d1, ...
	  rvecDeltaNorm_h2d2, rvecOmegaM1_h2d2, plot_marker_h2d2, 'color', plot_color_h2d2 );
	legend( ...
	  "H1D0", ...
	  "H1D1", ...
	  "H2D0", ...
	  "H2D1", ...
	  "H2D2", ...
	  "location", "north" );
	xlabel( "delta norm" );
	ylabel( "omega model 1" );
	title( "omega model 1 vs delta norm" );
	grid on;
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rvecDeltaNorm_h1d0, rvecOmegaM2_h1d0, plot_marker_h1d0, 'color', plot_color_h1d0, ...
	  rvecDeltaNorm_h1d1, rvecOmegaM2_h1d1, plot_marker_h1d1, 'color', plot_color_h1d1, ...
	  rvecDeltaNorm_h2d0, rvecOmegaM2_h2d0, plot_marker_h2d0, 'color', plot_color_h2d0, ...
	  rvecDeltaNorm_h2d1, rvecOmegaM2_h2d1, plot_marker_h2d1, 'color', plot_color_h2d1, ...
	  rvecDeltaNorm_h2d2, rvecOmegaM2_h2d2, plot_marker_h2d2, 'color', plot_color_h2d2 );
	legend( ...
	  "H1D0", ...
	  "H1D1", ...
	  "H2D0", ...
	  "H2D1", ...
	  "H2D2", ...
	  "location", "north" );
	xlabel( "delta norm" );
	ylabel( "omega model 2" );
	title( "omega model 2 vs delta norm" );
	grid on;
	%
	%
	%
	numColors = mygetfield( prm, "numColors", 1000 );
	sizeBigX = 201;
	sizeBigP = 203;
	bigXLo = bigX0 + min([ min(matDelta_h1d0(1,:)), min(matDelta_h1d1(1,:)) ]);
	bigXLo = min([ bigXLo, bigX_secret ]);
	bigXHi = bigX0 + max([ max(matDelta_h1d0(1,:)), max(matDelta_h1d1(1,:)) ]);
	bigXHi = max([ bigXHi, bigX_secret ]);
	bigPLo = bigP0 + min([ min(matDelta_h1d0(2,:)), min(matDelta_h1d1(2,:)) ]);
	bigPLo = min([ bigPLo, bigP_secret ]);
	bigPHi = bigP0 + max([ max(matDelta_h1d0(2,:)), max(matDelta_h1d1(2,:)) ]);
	bigPHi = max([ bigPHi, bigP_secret ]);
	rvecBigX = linspace( bigXLo-0.3*(bigXHi-bigXLo), bigXHi+0.3*(bigXHi-bigXLo), sizeBigX );
	rvecBigP = linspace( bigPLo-0.3*(bigPHi-bigPLo), bigPHi+0.3*(bigPHi-bigPLo), sizeBigP );
	[ matBigX, matBigP ] = meshgrid( rvecBigX, rvecBigP );
	dat_calcOmega = extFit_calcOmega_mat( matBigX, matBigP, rvecX, rvecF, rvecW );
	matOmegaAc = dat_calcOmega.matOmega;
	%
	rvecDeltaX = rvecBigX - bigX0;
	rvecDeltaP = rvecBigP - bigP0;
	[ matDeltaX, matDeltaP ] = meshgrid( rvecDeltaX, rvecDeltaP );
	matOmegaM1 = omega0 + matH1(2,1) * matDeltaX .* matDeltaP ....
	 + matDeltaX .* ( vecG(1) + 0.5*matH1(1,1)*matDeltaX ) ...
	 + matDeltaP .* ( vecG(2) + 0.5*matH1(2,2)*matDeltaP );
	matOmegaM2 = omega0 + matH2(2,1) * matDeltaX .* matDeltaP ....
	 + matDeltaX .* ( vecG(1) + 0.5*matH2(1,1)*matDeltaX ) ...
	 + matDeltaP .* ( vecG(2) + 0.5*matH2(2,2)*matDeltaP );
	%
	omegaMax = max(max(matOmegaAc));
	omegaMin = 0.0;
	%omegaMax = max(max(matOmegaM1));
	%omegaMin = min(min(matOmegaM1));
	funch_viz = @(om)( 1 + round((numColors-1)*cap( (om-omegaMin)/(omegaMax-omegaMin) ,0.0,1.0)) );
	%
	%
	%
	numFigs++; figure(numFigs);
	image( rvecBigX, rvecBigP, funch_viz(matOmegaAc) );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	plot( bigX_secret, bigP_secret, 'w*', 'markersize', 25, 'linewidth', 2 );
	plot( bigX0, bigP0, 'ws', 'markersize', 20, 'linewidth', 2 );
	plot( ...
	  bigX0+matDelta_h1d0(1,:), bigP0+matDelta_h1d0(2,:), plot_marker_h1d0, 'color', plot_color_h1d0, ...
	  bigX0+matDelta_h1d1(1,:), bigP0+matDelta_h1d1(2,:), plot_marker_h1d1, 'color', plot_color_h1d1, ...
	  bigX0+matDelta_h2d0(1,:), bigP0+matDelta_h2d0(2,:), plot_marker_h2d0, 'color', plot_color_h2d0, ...
	  bigX0+matDelta_h2d1(1,:), bigP0+matDelta_h2d1(2,:), plot_marker_h2d1, 'color', plot_color_h2d1, ...
	  bigX0+matDelta_h2d2(1,:), bigP0+matDelta_h2d2(2,:), plot_marker_h2d2, 'color', plot_color_h2d2 );
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omega vs bigX, bigP" );
	%
	numFigs++; figure(numFigs);
	image( rvecBigX, rvecBigP, funch_viz(matOmegaM1) );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	plot( bigX_secret, bigP_secret, 'w*', 'markersize', 25, 'linewidth', 2 );
	plot( bigX0, bigP0, 'ws', 'markersize', 20, 'linewidth', 2 );
	plot( ...
	  bigX0+matDelta_h1d0(1,:), bigP0+matDelta_h1d0(2,:), plot_marker_h1d0, 'color', plot_color_h1d0, ...
	  bigX0+matDelta_h1d1(1,:), bigP0+matDelta_h1d1(2,:), plot_marker_h1d1, 'color', plot_color_h1d1, ...
	  bigX0+matDelta_h2d0(1,:), bigP0+matDelta_h2d0(2,:), plot_marker_h2d0, 'color', plot_color_h2d0, ...
	  bigX0+matDelta_h2d1(1,:), bigP0+matDelta_h2d1(2,:), plot_marker_h2d1, 'color', plot_color_h2d1, ...
	  bigX0+matDelta_h2d2(1,:), bigP0+matDelta_h2d2(2,:), plot_marker_h2d2, 'color', plot_color_h2d2 );
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omega model 1 vs bigX, bigP" );
	%
	numFigs++; figure(numFigs);
	image( rvecBigX, rvecBigP, funch_viz(matOmegaM2) );
	axis("square");
	set( get(gcf,"children"), "ydir", "normal" );
	colormap(mycmap(numColors));
	hold on;
	plot( bigX_secret, bigP_secret, 'w*', 'markersize', 25, 'linewidth', 2 );
	plot( bigX0, bigP0, 'ws', 'markersize', 20, 'linewidth', 2 );
	plot( ...
	  bigX0+matDelta_h1d0(1,:), bigP0+matDelta_h1d0(2,:), plot_marker_h1d0, 'color', plot_color_h1d0, ...
	  bigX0+matDelta_h1d1(1,:), bigP0+matDelta_h1d1(2,:), plot_marker_h1d1, 'color', plot_color_h1d1, ...
	  bigX0+matDelta_h2d0(1,:), bigP0+matDelta_h2d0(2,:), plot_marker_h2d0, 'color', plot_color_h2d0, ...
	  bigX0+matDelta_h2d1(1,:), bigP0+matDelta_h2d1(2,:), plot_marker_h2d1, 'color', plot_color_h2d1, ...
	  bigX0+matDelta_h2d2(1,:), bigP0+matDelta_h2d2(2,:), plot_marker_h2d2, 'color', plot_color_h2d2 );
	hold off;
	grid on;
	xlabel( "bigX" );
	ylabel( "bigP" );
	title( "omega model 2 vs bigX, bigP" );
return;
%end
