function [ matHWB, datOut ] = getWellBehavedHessian( matH, prm=[] )
	datOut = [];
	datOut.muTot = 0.0;
	sz = size(matH,1);
	assert( isrealarray(matH,[sz,sz]) );
	assert( issymmetric(matH) );
	vecS = mygetfield( prm, "vecS", [] );
	if ( isempty(vecS) )
		vecS = ones(sz,1);
	endif
	assert( isrealarray(vecS,[sz,1]) );
	assert( 0.0 < min(vecS) );
	vecSInv = 1.0 ./ vecS;
	%
	matI = eye(sz,sz);
	if ( 0.0 == max(max(abs(matH))) )
		msg( __FILE__, __LINE__, "Warning: Input Hessian is zero." );
		matH = matI;
	endif
	matHS = (matH.*vecSInv).*(vecSInv'); % Autobroadcast.
	clear matH; % To avoid accidental use. Move this if you need matH later.
	%
	% Make sure Hessian is positive-definite.
	[ matR, cholFlag ] = chol( matHS );
	if ( 0 ~= cholFlag || min(diag(matR)) <= sqrt(eps)*max(abs(diag(matR))) )
		% Matrix is not clearly positive-definite.
		vecEig = eig( matHS );
		muCrit = max([ -min(vecEig), 0.0 ]);
		muSafe = muCrit + sqrt(eps)*max(max(abs(matHS)));
	else
		muSafe = 0.0;
	endif
	matHSPD = matHS + muSafe*matI;
	datOut.muTot += muSafe;
	%
	% Were we asked to ensure the minimum value is not below some allowed value?
	f0 = mygetfield( prm, "f0", [] );
	if ( isempty(f0) )
		% ... No, we were not.
		matHWB = (matHSPD.*vecS).*(vecS');
		return;
	endif
	% ... Yes we were.
	assert( isrealscalar(f0) );
	vecG0 = mygetfield( prm, "vecG0", [] );
	if ( isempty(vecG0) )
		msg( __FILE__, __LINE__, "ERROR: Handling for fMinAllowed requires both f0 and vecG0." );
		error( "Handling for fMinAllowed requires both f0 and vecG0." );
	endif
	assert( isrealarray(vecG0,[sz,1]) );
	vecGS = vecG0.*vecSInv;
	clear vecG0;
	fMinAllowed = mygetfield( prm, "fMinAllowed", -0.1*f0 );
	assert( isrealscalar(fMinAllowed) );
	assert( fMinAllowed < f0 );
	%
	% We need to find some mu such that f0 - 0.5*vecG'*inv(matHSPD+mu*matI)*vecG >= fMin.
	% This corresponds to f = f0 + vecDelta'*vecGS + 0.5*vecDelta'*(matHSPD+mu*matI)*vecDelta
	%  with vecDelta = -( matHSPD + mu*matI ) \ vecG.
	% Note that vecG'*inv(matA)*vecG = sumsq(chol(matA)'\vecG)^2... Not that such is used in this version.
	funchFNewt = @(mu)( funcFNewt( mu, matHSPD, f0, vecGS, matHS ) );
	if ( funchFNewt( 0.0 ) >= fMinAllowed )
		% No need to make any further modification.
		matHWB = (matHSPD.*vecS).*(vecS');
		return;
	endif
	% This muHi is just speculative.
	%msg( __FILE__, __LINE__, "vvv" );
	hLo = sqrt(sum(sum(matHS.^2,1)))/(sz*sz);
	hHi = max(max(abs(matHS)));
	gLo = sqrt(sum(vecGS.^2))/sz;
	gHi = max(max(abs(vecGS)));
	fLo = sqrt( f0^2 );
	fHi = abs(f0) + abs(fMinAllowed) + hHi*(gHi^2);
	%deltaNormHi = gHi / hLo;
	muHi = 10.0*fHi*hHi/(gLo^2);
	%bigDelta = norm(vecGS) / sqrt(max(sum(matHS.^2,1)));
	%bigF = abs(f0) + abs(fMinAllowed) + sz*max(abs(vecGS));
	%muHi_old = 100.0*(abs(f0)+abs(fMinAllowed))*max(max(abs(matHS)))/sumsq(vecGS)
	%muHi = 1E5
	muBracket = [ 0.0, muHi ];
	ffnBracket = [ funchFNewt(0.0), funchFNewt(muHi) ];
	% fzero() does not utilize analytic derivatives, which are redaily available in this case.
	%  d/ds (M^-1) = - M^-1 * (d/ds M) * M^-1.
	muStayPositive = fzero( funchFNewt, [ 0.0, muHi ] );
	datOut.muTot += muStayPositive;
	matHWB = (( matHSPD + muStayPositive*matI ).*vecS).*(vecS');
	%msg( __FILE__, __LINE__, "^^^" );
return;
endfunction

function f = funcFNewt( mu, matHSPD, f0, vecGS, matHS )
	matR = chol( matHSPD + mu*eye(size(matHSPD)) );
	vecDelta = matR \ ( matR' \ (-vecGS) );
	f = f0 + vecDelta'*vecGS + 0.5*(vecDelta'*matHS*vecDelta);
return;
endfunction
