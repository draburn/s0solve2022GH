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
		%setprngstates(38566912);
		setprngstates(69539920);
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
	[ muLim, retCode, datOut ] = extFit_findMuLim( omegaLim, omega0, vecG, matH, matR, prm );
	%
	rvecMu = datOut.muCrit + datOut.muScale * 100.0 * linspace(0.0,1.0,101).^4;
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
	grid on;
	%
	rvecMask = ( rvecOmega > omegaLim - abs( omega0 - omegaLim ) );
	numFigs++; figure(numFigs);
	plot( ...
	  rvecMu(rvecMask), rvecOmega(rvecMask), 'o-', ...
	  [min(rvecMu), max(rvecMu)], omegaLim*[1,1], 'k-', ...
	  muLim, omegaLim, '*', 'markersize', 30 );
	grid on
	%
	%numFigs++; figure(numFigs);
	%plot( ...
	%  rvecMu, rvecDeltaNorm, 'o-' );
	%grid on;
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rvecMu(rvecMask), rvecDeltaNorm(rvecMask), 'o-' );
	grid on;
	%
	%numFigs++; figure(numFigs);
	%plot( ...
	%  rvecDeltaNorm, rvecOmega, 'o-' );
	%grid on;
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rvecDeltaNorm(rvecMask), rvecOmega(rvecMask), 'o-' );
	grid on;

