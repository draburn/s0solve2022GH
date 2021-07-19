% function [ xOfCand, meritOfCand, datOut ] = extFinder( xVals, fVals, prm = [], datIn = [] );
	setprngstates(0);
	xVals = sort(randn(1,200));
	fVals = abs(xVals).^4
	prm = [];
	datIn = [];
	%
	%plot( xVals, fVals, 'o-' );
	%grid on;
	%
	extFinder__init;
	thisFile = "extFinder";
	if (haveCand)
		thisFile = [ "RETURN from " thisFile ];
		return;
	end
	%
	%
	bigS0 = mygetfield( prm, "bigS0", [] );
	bigP0 = mygetfield( prm, "bigP0", [] );
	if ( isempty(bigS0) || isempty(bigP0) )
		prm_dlogfdf = mygetfield( prm, "prm_dlogfdf", [] );
		[ bigX_temp, bigP_temp ] = extFinder_dlogdf( xVals, fVals, prm_dlogfdf );
		if ( isempty(bigP0) )
			bigP0 = bigP_temp;
		end
		assert( isrealscalar(bigP0) );
		if ( isempty(bigS0) )
			bigP = bigP0;
			bigS = bigX_temp;
			[ bigA, bigB, bigC ] = extFinder_getFit( bigS, bigP, xVals, fVals, nOfPtWiseMin );
			xExt = bigS - bigDelta*mypow( bigB/(bigA*bigP), 1.0/(bigP-1.0) )
			assert( xVals(nOfPtWiseMin-1) < xExt ); % What to do if this fails?!?!
			assert( xExt < xVals(nOfPtWiseMin+1) ); % What to do if this fails?!?!
			bigS0 = bigX_temp;
		end
		assert( isrealscalar(bigS0) );
		clear bigX_temp;
		clear bigP_temp;
		clear prm_dlogfdf;
		%
		%
		echo__bigS0 = bigS0
		echo__bigP0 = bigP0
	end
	assert( isrealscalar(bigS0) );
	assert( isrealscalar(bigP0) );
	bigS = bigS0;
	bigP = bigP0;
	[ bigA, bigB, bigC ] = extFinder_getFit( bigS, bigP, xVals, fVals, nOfPtWiseMin );
	%
	%
	funchGModel = @(x)( bigG0 + bigG1*( ...
	  bigA*abs((x-bigS)/bigDelta).^bigP + bigB*(x-bigS)/bigDelta + bigC ) );
	plot( ...
	  xVals, gVals, 'o-', ...
	  xVals, funchGModel(xVals), 'x-' );
	grid on;
	%
	tic();
	extFinder_getOmega( bigS, bigP, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( bigS, 7.5, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( 0.0, bigP, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( 0.0, 7.5, xVals, fVals, nOfPtWiseMin )
	%
	epsS = 1e-5;
	epsP = 1e-2;
	omegaMM = extFinder_getOmega( bigS-epsS, bigP-epsP, xVals, fVals, nOfPtWiseMin );
	omegaM0 = extFinder_getOmega( bigS-epsS, bigP,      xVals, fVals, nOfPtWiseMin );
	omegaMP = extFinder_getOmega( bigS-epsS, bigP+epsP, xVals, fVals, nOfPtWiseMin );
	omega0M = extFinder_getOmega( bigS,      bigP-epsP, xVals, fVals, nOfPtWiseMin );
	omega00 = extFinder_getOmega( bigS,      bigP,      xVals, fVals, nOfPtWiseMin );
	omega0P = extFinder_getOmega( bigS,      bigP+epsP, xVals, fVals, nOfPtWiseMin );
	omegaPM = extFinder_getOmega( bigS+epsS, bigP-epsP, xVals, fVals, nOfPtWiseMin );
	omegaP0 = extFinder_getOmega( bigS+epsS, bigP,      xVals, fVals, nOfPtWiseMin );
	omegaPP = extFinder_getOmega( bigS+epsS, bigP+epsP, xVals, fVals, nOfPtWiseMin );
	%
	omegaDS = ( omegaP0 - omegaM0 ) / (2.0*epsS)
	omegaDP = ( omega0P - omega0M ) / (2.0*epsP)
	omegaDDS = ( omegaP0 + omegaM0 - 2*omega00 ) / ( epsS^2 )
	omegaDDP = ( omega0P + omega0M - 2*omega00 ) / ( epsP^2 )
	%
	deltaS = -omegaDS/omegaDDS
	deltaP = -omegaDP/omegaDDP
	extFinder_getOmega( bigS + 1.0*deltaS, bigP + 1.0*deltaP, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( bigS + 0.9*deltaS, bigP + 0.9*deltaP, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( bigS + 0.8*deltaS, bigP + 0.8*deltaP, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( bigS + 0.5*deltaS, bigP + 0.5*deltaP, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( bigS + 0.2*deltaS, bigP + 0.2*deltaP, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( bigS + 0.1*deltaS, bigP + 0.1*deltaP, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( bigS + 0.01*deltaS, bigP + 0.01*deltaP, xVals, fVals, nOfPtWiseMin )
	extFinder_getOmega( bigS, bigP, xVals, fVals, nOfPtWiseMin )
	toc();
	bigS += 0.2*deltaS
	bigP += 0.2*deltaP
	%
	%
	extFinder_viz( xVals, fVals, -0.02, 0.02, 3.8, 4.2 );
return;
