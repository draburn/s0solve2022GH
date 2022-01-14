	clear;
	thisFile = "test_extFit_findMuLim"
	numFigs = 0;
	if (0)
		omegaLim = 0.5;
		omega0 = 1.0;
		vecG = [1;1];
		matH = [-1,0;0,-1];
		matR = [1,0;0,1];
	else
		setprngstates();
		%setprngstates(44975936);
		%setprngstates(38883696);
		%setprngstates(84427936); % Non-monotonic. It is known.
		%setprngstates(16576176);
		%setprngstates(82450176);
		%setprngstates(40131696);
		%setprngstates(38566912);
		%setprngstates(69539920);
		%setprngstates(82122480);
		%setprngstates(53282288);
		omegaLim = 50.0;
		omega0 = 100.0;
		vecG = randn(2,1);
		matH = randn(2,2);
		matH = matH' + matH;
		matR = randn(2,2);
		matR = matR'*matR;
	end
	prm = [];
	%
	mu = extFit_findMuOfOmega( omegaLim, omega0, vecG, matH, matR, prm );
	return;
	[ muLim, retCode, datOut ] = extFit_findMuLim( omegaLim, omega0, vecG, matH, matR, prm );
	vecDeltaLim = -(matH + muLim*matR)\vecG;
	%
	vecDelta = -(matH + muLim*matR) \ vecG;
	omegaAtMuLim = omega0 + (vecDelta'*vecG) + (0.5*vecDelta'*matH*vecDelta);
	rvecMu = 0.5*(datOut.muCrit+muLim) + datOut.muScale * 1000.0 * linspace(0.0,1.0,101).^4;
	for n=1:max(size(rvecMu))
		mu = rvecMu(n);
		vecDelta = -(matH + mu*matR) \ vecG;
		omega = omega0 + (vecDelta'*vecG) + (0.5*vecDelta'*matH*vecDelta);
		rvecDeltaNorm(n) = norm(vecDelta);
		rvecOmega(n) = omega;
	end
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rvecMu, rvecOmega, 'o-', ...
	  [min(rvecMu), max(rvecMu)], omegaLim*[1,1], 'k-', ...
	  muLim, omegaLim, '*', 'markersize', 30 );
	xlabel( "mu" );
	ylabel( "omega" );
	title( "omega vs mu (unmasked)" );
	grid on;
	%
	%rvecMask = ( rvecOmega > omegaLim - abs( omega0 - omegaLim ) );
	rvecMask = ( ones(size(rvecOmega)) > zeros(size(rvecOmega)) );
	numFigs++; figure(numFigs);
	plot( ...
	  rvecMu(rvecMask), rvecOmega(rvecMask), 'o-', ...
	  [min(rvecMu), max(rvecMu)], omegaLim*[1,1], 'k-', ...
	  muLim, omegaLim, '*', 'markersize', 30 );
	xlabel( "mu" );
	ylabel( "omega" );
	title( "omega vs mu (masked)" );
	grid on
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rvecMu(rvecMask), rvecDeltaNorm(rvecMask), 'o-' );
	grid on;
	xlabel( "mu" );
	ylabel( "||vecDelta||" );
	title( "||vecDelta|| vs mu (masked)" );
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rvecDeltaNorm(rvecMask), rvecOmega(rvecMask), 'o-', ...
	  [0.0, max(rvecDeltaNorm(rvecMask))], omegaLim*[1.0,1.0], 'k-', ...
	  norm(vecDeltaLim), omegaAtMuLim, '*', 'markersize', 30 );
	grid on;
	xlabel( "||vecDelta||" );
	ylabel( "omega" );
	title( "omega vs ||vecDelta|| (masked)" );
