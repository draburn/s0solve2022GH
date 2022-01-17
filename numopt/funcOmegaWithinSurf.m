% Function...


function [ omegaVals, vecNablaOmegaVals ] = funcOmegaWithinSurf( vecXVals, funchOmegaBase, funchSurf, prm=[] )
	if ( 3 > nargin || 4 < nargin )
		msg( __FILE__, __LINE__, "Bad nargin." );
		print_usage();
		return; % Superfluous?
	elseif ( 2 < nargout )
		msg( __FILE__, __LINE__, "Bad nargout." );
		print_usage();
		return; % Superfluous?
	end
	%
	h0 = mygetfield( prm, "h0", 0.01 );
	tau = mygetfield( prm, "tau", 0.01 );
	numVals = size(vecXVals,2);
	%
	if ( 1 < numVals )
		[ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = funchSurf( vecXVals );
		vecDVals = vecXVals - vecSVals;
		switch (nargout)
		case 1
			omegaVals_i = funchOmegaBase( vecXVals );
			[ omegaVals_s, vecNablaOmegaVals_s ] = funchOmegaBase( vecSVals );
			omegaVals_o = omegaVals_s + sum( vecDVals .* vecNablaOmegaVals_s, 1 ) ...
			  + 0.5 * sumsq( vecDVals, 1 ) .* ( h0 + sqrt(sumsq( vecNablaOmegaVals_s, 1 )) / tau );
			outFlagVals = ( sum(vecDVals.*vecNHatVals,1) > 0.0 );
			omegaVals = omegaVals_i;
			omegaVals(outFlagVals) = omegaVals_o(outFlagVals);
		case 2
			[ omegaVals_i, vecNablaOmegaVals_i, matNabla2OmegaVals_i ] = funchOmegaBase( vecXVals );
			[ omegaVals_s, vecNablaOmegaVals_s, matNabla2OmegaVals_s ] = funchOmegaBase( vecSVals );
			omegaVals_o = omegaVals_s + sum( vecDVals .* vecNablaOmegaVals_s, 1 ) ...
			  + 0.5 * sumsq( vecDVals, 1 ) .* ( h0 + sqrt(sumsq( vecNablaOmegaVals_s, 1 )) / tau );
			outFlagVals = ( sum(vecDVals.*vecNHatVals,1) > 0.0 );
			omegaVals = omegaVals_i;
			omegaVals(outFlagVals) = omegaVals_o(outFlagVals);
			msg( __FILE__, __LINE__, "The rest of this case is not implemented." );
			error ( "To-do." );
		otherwise
		error( "Impossible case." );
		end
	return;
	end
	%
	% numVals is 1.
	[ vecS, vecNHat, vecUHat, matNablaST ] = funchSurf( vecXVals );
	vecD = vecXVals - vecS;
	if ( vecD'*vecNHat <= 0.0 )
		% We're inside.
		switch (nargout)
		case 1
			omegaVals = funchOmegaBase( vecXVals );
		case 2
			[ omegaVals, vecNablaOmegaVals ] = funchOmegaBase( vecXVals );
		otherwise
			error( "Impossible case." );
		end
	return;
	end
	%
	% numVals is 1 and we're outside.
	switch (nargout)
	case 1
		[ omega_s, vecNablaOmega_s ] = funchOmegaBase( vecS );
		omegaVals = omega_s + vecD'*vecNablaOmega_s + 0.5 * sumsq(vecD) * ( h0 + norm(vecNablaOmega_s)/tau );
	case 2
		[ omega_s, vecNablaOmega_s, matNabla2Omega_s ] = funchOmegaBase( vecS );
		omegaVals = omega_s + svecD'*vecNablaOmega_s + 0.5 * sumsq(vecD) * ( h0 + norm(vecNablaOmega_s)/tau );
		msg( __FILE__, __LINE__, "The rest of this case is not implemented." );
		error ( "To-do." );
	otherwise
		error( "Impossible case." );
	end
return;
end


%!test
%!	msg( __FILE__, __LINE__, "Performing basic execution test." );
%!	%msg( __FILE__, __LINE__, "SKIPPING TEST." ); return;
%!	setprngstates(0);
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = 2 + round(2.0*abs(randn()));
%!	numVals = 50 + round(10.0*abs(randn()));
%!	%
%!	vecXCent_surf = randn(sizeX,1);
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent_surf );
%!	vecXCent_base = randn(sizeX,1);
%!	funchOmegaBase = @(dummyX) funcOmegaEllip( dummyX, vecXCent_base );
%!	%
%!	vecXVals = vecXCent_surf + 0.4*randn(sizeX,numVals);
%!	%
%!	%
%!	% Test with minimal input.
%!	for n=1:numVals
%!		omega = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf );
%!		assert( isrealscalar(omega) );
%!		%[ omega, vecNablaOmega ] = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf );
%!		%assert( isrealscalar(omega) );
%!		%assert( isrealarray(vecNablaOmega,[sizeX,1]) );
%!	end
%!	%
%!	%
%!	% Test with maximal input.
%!	bigR = 0.01 + abs(randn);
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA_surf = (0.01*eye(sizeX,sizeX)) + (matA0'*matA0);
%!	debugMode_surf = true;
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent_surf, bigR, matA_surf, debugMode_surf );
%!	%
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA_base = (0.01*eye(sizeX,sizeX)) + (matA0'*matA0);
%!	omega0 = abs(randn);
%!	omega1 = 0.01 + abs(randn);
%!	debugMode_base = true;
%!	funchOmega = @(dummyX) funcOmegaEllip( dummyX, vecXCent_base, matA_base, omega0, omega1, debugMode_base );
%!	%
%!	prm = [];
%!	omegaVals = funcOmegaWithinSurf( vecXVals, funchOmegaBase, funchSurf, prm );
%!	assert( isrealarray(omegaVals,[1,numVals]) );
%!	%[ omegaVals, vecNablaOmegaVals ] = funcOmegaWithinSurf( vecXVals, funchOmegaBase, funchSurf, prm );
%!	%assert( isrealarray(omegaVals,[1,numVals]) );
%!	%assert( isrealarray(vecNablaOmegaVals,[sizeX,numVals]) );
%!	%
%!	%
%!	% Test vectorization.
%!	epsOmega = sqrt(eps*sumsq(reshape(omegaVals,[],1)))/numVals;
%!	%epsNablaOmega = sqrt(eps*sumsq(reshape(vecNablaOmegaVals,[],1)))/numVals;
%!	for n=1:numVals
%!		omega = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf, prm );
%!		assert( reldiff(omega,omegaVals(n),epsOmega) < sqrt(eps) );
%!		%[ omega, vecNablaOmega ] = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf );
%!		%assert( reldiff(omega,omegaVals(n),epsOmega) < sqrt(eps) );
%!		%assert( reldiff(vecNablaOmega,vecNablaOmegaVals(:,n),epsNablaOmega) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );