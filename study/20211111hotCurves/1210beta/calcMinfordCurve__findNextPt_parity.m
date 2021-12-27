function [ vecX, datOut ] = calcMinfordCurve__findNextPt( funchOmega, funchG, onSurf0, vecX0, vecXC, bigR, matS=[], prm=[] )
	thisFile = "calcMinfordCurve__findNextPt";
	%msg( thisFile, __LINE__, "THIS IS A WORK-IN-PROGRESS!" );
	%
	if (~isempty(onSurf0))
		assert(isscalar(onSurf0));
	end
	sizeX = size(vecX0,1);
	assert(isrealarray(vecX0,[sizeX,1]));
	assert(isrealarray(vecXC,[sizeX,1]));
	assert(isrealscalar(bigR));
	assert( 0.0 < bigR );
	if (~isempty(matS))
		assert(isrealarray(matS,[sizeX,sizeX]));
		matS_nonEmpty = matS;
	else
		matS_nonEmpty = eye(sizeX,sizeX);
		% Use this for convenience, but, for efficiency,
		% let sub-modules know that matS is actually empty.
	end
	datOut = [];
	%
	% Validate onSurf0 wrt vecX, etc.
	% vecS = matS*(vecX-vecXC).
	% If s < R, we must be inside the surface.
	% If s ~ R, we must rely on onSurf.
	% If s > R, something has gone wrong.
	s0 = norm( matS_nonEmpty*(vecX0-vecXC) );
	if (isempty(onSurf0))
		if ( s0 < bigR )
			onSurf0 = false;
		else
			onSurf0 = true;
			vecX0 = vecXC + (vecX0-vecXC)*(bigR/s0);
		end
	elseif (onSurf0)
		assert( abs( s0 - bigR ) <= sqrt(eps)*bigR );
	else
		assert( s0 < bigR );
	end
	%
	%
	%
	omega0 = funchOmega( vecX0 );
		foo = funchG( vecX0 ); %%% FORCE "PARITY".
	assert(isrealscalar(omega0));
	assert( omega0 >= 0.0 );
	%
	%
	useHighAccuracyParams = false;
	if (useHighAccuracyParams)
		iterLimit = 5000;
		gActualAbsNormTol = 1e-12;
		gSurfAbsNormTol = 1e-10;
		gSurfRelNormTol = 1e-4;
		deltaNormTarget = 1e-2;
		btLimit = 50;
		btFactor = 0.7;
		deltaNormTol = 1e-4;
	else
		iterLimit = 100;
		gActualAbsNormTol = 1e-8;
		gSurfAbsNormTol = 1e-7;
		gSurfRelNormTol = 1e-3;
		deltaNormTarget = 1e-1;
		btLimit = 5;
		btFactor = 0.1;
		deltaNormTol = 1e-3;
	end
	%
	%
	%
	iterCount = 0;
	onSurf = onSurf0;
	vecX = vecX0;
	omega = omega0;
	while (1)
		vecGActual = funchG( vecX );
			foo = funchOmega( vecX ); %%% FORCE "PARITY".
		assert(isrealarray(vecGActual,[sizeX,1]));
		s = norm( matS_nonEmpty*(vecX-vecXC) );
		if (onSurf)
			assert( abs( s - bigR ) <= sqrt(eps)*bigR );
			vecD = vecX-vecXC;
			vecGBall = (bigR/s)*( vecGActual - (vecD'*vecGActual/s^2)*(matS_nonEmpty'*matS_nonEmpty*vecD) );
			assert( isrealarray(vecGBall,[sizeX,1]) );
			if ( norm(vecGBall) <= gSurfAbsNormTol )
				%msg( thisFile, __LINE__, "Reached gSurfAbsNormTol." );
				return;
			end
			if ( norm(vecGBall) <= gSurfRelNormTol*norm(vecGActual) )
				%msg( thisFile, __LINE__, "Reached gSurfRelNormTol." );
				return;
			end
		else
			assert( s < bigR );
			vecGBall = vecGActual;
			if ( norm(vecGActual) <= gActualAbsNormTol )
				%msg( thisFile, __LINE__, "Reached gActualAbsNormTol (inside surface)." );
				return;
			end
		end
		%
		iterCount++;
		if ( iterCount > iterLimit )
			msg( thisFile, __LINE__, "Reached iterLimit." );
			msg( thisFile, __LINE__, "(Consider implementing a more sophisticated solver to be used conditionally.)" );
			return;
		end
		%
		%
		% Fall off surface if gradient-descent direction is inward.
		if (onSurf)
		if ( vecGActual'*(vecX-vecXC) > 0.0 )
			%msg( thisFile, __LINE__, "Fell off surface!" );
			if (0)
				msg( thisFile, __LINE__, "" );
				msg( thisFile, __LINE__, "vvvvvvvvvv" );
				echo__vecXC = vecXC
				echo__bigR = bigR
				echo__matS = matS
				echo__vecX = vecX
				msg( thisFile, __LINE__, "^^^^^^^^^^" );
				msg( thisFile, __LINE__, "" );
			end
			%
			onSurf = false;
			if ( norm(vecGActual) <= gActualAbsNormTol )
				msg( thisFile, __LINE__, "Reached gActualAbsNormTol (upon falling off surface)." );
				return;
			end
			vecGBall = vecGActual;
		end
		end
		%
		%
		%
		if (onSurf)
			%deltaT = deltaNormTarget/norm(vecGBall); % Old
			%
			%
			vecV = -vecGBall;
			epsFD = (0.1*deltaNormTol)/norm(vecV);
			%
			% vecXOff would represent a vecX that is off the surface.
			% vecDOff = vecXOff - vecXC.
			vecDOffP = vecX + (vecV*epsFD) - vecXC;
			vecXP = vecXC + (vecDOffP*(bigR/norm(matS_nonEmpty*vecDOffP)));
			omegaP = funchOmega(vecXP);
				foo = funchG( vecXP ); %%% FORCE "PARITY".
			assert( isrealscalar(omegaP) );
			assert( 0.0 <= omegaP );
			omega0 = omega;
			dOmega = vecGBall'*vecV;
			% omegaModel = omega0 + dOmega*t + c*(t^2);
			% d2Omega = 2.0*c;
			d2Omega = 2.0*( omegaP - omega0 - (dOmega*epsFD) ) / (epsFD^2);
			%
			deltaT = calcLinishRootOfQuad( 0.5*d2Omega, dOmega, omega0 );
			assert( 0.0 < deltaT );
			%
			if (0)
				numTVals = 1001;
				tVals = linspace(-1.0,1.0,numTVals);
				omegaVals = zeros(1,numTVals);
				for n=1:numTVals
					vecDOff = vecX + (vecV*tVals(n)) - vecXC;
					omegaVals(n) = funchOmega( vecXC + (vecDOff*(bigR/norm(matS_nonEmpty*vecDOff))) );
						foo = funchG( vecXC + (vecDOff*(bigR/norm(matS_nonEmpty*vecDOff))) ); %%% FORCE "PARITY".
				end
				plot( ...
				  tVals, omegaVals, 'o-', ...
				  tVals, omega0+(dOmega*tVals)+(0.5*d2Omega*(tVals.^2)), 'x-', ...
				  deltaT+(0*tVals), omegaVals, '^-', ...
				  0*tVals, omegaVals, 'kv-' );
				grid on;
				assert(0);
			end
			%
			%
			onSurf_next = true; % Until and unless fall off at start of next iter.
			vecX_next = vecX - (vecGBall*deltaT);
			vecX_next = vecXC + (vecX_next-vecXC)*(bigR/norm(matS_nonEmpty*(vecX_next-vecXC)));
			omega_next = funchOmega( vecX_next );
				foo = funchG( vecX_next ); %%% FORCE "PARITY".
			assert( isrealscalar(omega_next) );
			assert( 0.0 <= omega_next );
			%
			btCount = 0;
			while ( omega_next >= omega )
				btCount++;
				if ( btCount >= btLimit )
					%msg( thisFile, __LINE__, "Reached backtracking limit (while on surface)." );
					% To force the caller to stop, consider: vecX = vecX0;
					return;
				end
				if ( norm(vecX_next-vecX) <= deltaNormTol )
					%msg( thisFile, __LINE__, "Reached deltaNormTol (while on surface)." );
					return;
				end
				%
				%msg( thisFile, __LINE__, sprintf( "Backtracking %d...", btCount ) );
				deltaT *= btFactor;
				vecX_next = vecX - (vecGBall*deltaT);
				vecX_next = vecXC + (vecX_next-vecXC)*(bigR/norm(matS_nonEmpty*(vecX_next-vecXC)));
				omega_next = funchOmega( vecX_next );
					foo = funchG( vecX_next ); %%% FORCE "PARITY".
				assert( isrealscalar(omega_next) );
				assert( 0.0 <= omega_next );
			end
		else
			%deltaT = deltaNormTarget/norm(vecGBall); % OLD
			%
			%
			vecV = -vecGBall;
			epsFD = (0.1*deltaNormTol)/norm(vecV);
			%
			vecXP = vecX + (vecV*epsFD);
			omegaP = funchOmega(vecXP);
				foo = funchG(vecXP); %%% FORCE "PARITY".
			assert( isrealscalar(omegaP) );
			assert( 0.0 <= omegaP );
			omega0 = omega;
			dOmega = vecGBall'*vecV;
			% omegaModel = omega0 + dOmega*t + c*(t^2);
			c = ( omegaP - omega0 - (dOmega*epsFD) ) / (epsFD^2);
			%
			deltaT = calcLinishRootOfQuad( c, dOmega, omega0 );
			assert( 0.0 < deltaT );
			%
			%
			onSurf_next = false; % Unless changed below.
			vecX_next = vecX - (vecGBall*deltaT);
			%
			% Check if went past surface.
			if ( norm(matS_nonEmpty*(vecX_next-vecXC)) >= bigR )
				% Reduce deltaT to hit surface.
				%msg( thisFile, __LINE__, "Hit surface." );
				onSurf_next = true; % Unless changed below, because the on surface point is bad.
				vecSXMXC = matS_nonEmpty*(vecX-vecXC);
				vecSG = matS_nonEmpty*vecGBall;
				c2 = vecSG'*vecSG;
				c1 = -2.0*(vecSG'*vecSXMXC);
				c0 = vecSXMXC'*vecSXMXC - bigR^2;
				%
				% getQuadGoodPt() here is wrong;
				% we always want the "forward" solution,
				% and we have no guarantee about the sign of c1.
				%deltaT = getQuadGoodPt( c2, c1, c0 );
				% However, the problem is simple
				% since there should be exactly one positive solution...
				% If we are already at the solution then deltaT = 0.
				assert( c2 > 0.0 );
				assert( c0 < (eps^0.75)*(bigR^2) );
				d = (c1^2) - (4.0*c0*c2);
				assert( d > -(eps^0.75)*( (c1^2) + abs(c0*c2) ) );
				if ( d > 0.0 )
					sqrtD = sqrt(d);
				else
					sqrtD = 0.0;
				end
				deltaT = ( -c1 + sqrtD ) / ( 2.0*c2 );
				assert( deltaT >= -(eps^0.75)*(abs(c1/c2)+abs(sqrtD/c2)) );
				if ( 0.0 > deltaT )
					deltaT = 0.0;
				end
				%
				vecX_next = vecX - (vecGBall*deltaT);
			end
			%
			omega_next = funchOmega( vecX_next );
				foo = funchG( vecX_next ); %%% FORCE "PARITY".
			assert( isrealscalar(omega_next) );
			assert( 0.0 <= omega_next );
			%
			btCount = 0;
			while ( omega_next >= omega )
				btCount++;
				if ( btCount >= btLimit )
					%msg( thisFile, __LINE__, "Reached backtracking limit (inside surface)." );
					return;
				end
				if ( norm(vecX_next-vecX) <= deltaNormTol )
					%msg( thisFile, __LINE__, "Reached deltaNormTol (inside surface)." );
					return;
				end
				%msg( thisFile, __LINE__, sprintf( "Backtracking %d...", btCount ) );
				%
				deltaT *= btFactor;
				onSurf_next = false; % We are definitely inside the surface now.
				vecX_next = vecX - (vecGBall*deltaT);
				omega_next = funchOmega( vecX_next );
					foo = funchG( vecX_next ); %%% FORCE "PARITY".
				assert( isrealscalar(omega_next) );
				assert( 0.0 <= omega_next );
			end
		end
		assert( omega_next < omega );
		onSurf_prev = onSurf;
		%
		vecX = vecX_next;
		omega = omega_next;
		onSurf = onSurf_next;
	end
end
