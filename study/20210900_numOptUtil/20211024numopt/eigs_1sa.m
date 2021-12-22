function [ lambdaSA, datOut ] = eigs_1sa( matH, vecG=[], prm = [] )
	thisFile = "eigs_1sa";
	eps050 = sqrt(eps);
	eps025 = sqrt(eps050);
	datOut = [];
	%
	probSize = size(matH,1);
	assert( probSize >= 1 );
	if (isempty(vecG) )
		vecG = ones(probSize,1);
	end
	assert( isrealarray(vecG,[probSize,1]) );
	assert( isrealarray(matH,[probSize,probSize]) );
	assert( issymmetric(matH) );
	%
	hScale = eigs( matH, 1 );
	matI = eye( probSize, probSize );
	%
	muFail = mygetfield( prm, "muFail", 0.0 );
	muSafe = mygetfield( prm, "muSafe", hScale*(1.0+eps025) );
	assert( isrealscalar(muFail) );
	assert( isrealscalar(muSafe) );
	assert( 0.0 <= muFail );
	assert( muFail < muSafe );
	%
	mu = muSafe;
	matM = matH + (mu*matI);
	[ matR, cholFlag ] = chol( matM );
	if ( 0 ~= cholFlag )
		error( "chol() failed for muSafe." );
	end
	clear opts;
	opts.tol = 0.1;
	opts.maxit = 10;
	opts.v0 = vecG;
	opts.cholB = true;
	[ vecPsi, x, eigsFlag ] = eigs( matI, matR, 1, 'lm', opts );
	if ( 0 ~= eigsFlag )
		error( "eigs() failed for muSafe" );
	end
	muSafe = mu;
	muCrit = mu - (1.0./x);
	muNext = muCrit*(1.0+eps025);
	if ( muNext <= muFail )
		msg( thisFile, __LINE__, "" );
		muNext = (muFail+muSafe)/2.0
	elseif ( muNext >= muSafe )
		msg( thisFile, __LINE__, "" );
		muNext = (muFail+muSafe)/2.0
	end
	%
	%
	%
	iterLimit = mygetfield( prm, "iterLimit", 10 );
	iterCount = 0;
	while (1)
		muPrev = mu;
		mu = muNext;
		iterCount++;
		msg( thisFile, __LINE__, sprintf( ...
		  "  %3d;  %10.3e, %10.3e;  %10.3e;  %10.3e, %10.3e.", ...
		  iterCount, muFail, muSafe, muCrit, muPrev, muNext ) );
		if ( iterCount > iterLimit )
			break;
		end
		assert( mu > muFail );
		assert( mu < muSafe );
		muLo = (1.0-eps025)*muFail + eps025*muSafe
		muHi = eps025*muFail + (1.0-eps025)*muSafe
		%
		matM = matH + (mu*matI);
		[ matR, cholFlag ] = chol( matM );
		if ( 0 == cholFlag )
			clear opts;
			opts.tol = 0.1;
			opts.maxit = 10;
			opts.v0 = vecG;
			opts.cholB = true;
			[ vecPsi, x, eigsFlag ] = eigs( matI, matR, 1, 'lm', opts );
			if ( 0 == eigsFlag )
				msg( thisFile, __LINE__, "" );
				muSafe = mu;
				muCrit = mu - (1.0./x)
				if ( abs(muCrit-muSafe) < eps050*hScale )
					msg( thisFile, __LINE__, "" );
					muNext = -1.0;
					break;
				end
				muNext = muCrit*( 1.0 + 2.0*opts.tol )
				if ( muNext < muFail*(1.0-eps025) )
					msg( thisFile, __LINE__, "" );
					muNext = 0.9*muFail+0.1*muSafe;
				elseif ( muNext > muSafe*(1.0+eps025) )
					msg( thisFile, __LINE__, "" );
					muNext = 0.9*muFail+0.1*muSafe;
				else
					msg( thisFile, __LINE__, "" );
					muNext = cap( muNext, muLo, muHi );
				end
			else
				msg( thisFile, __LINE__, "" );
				muFail = mu;
				muNext = (muFail+muSafe)/2.0;
			end
		else
			msg( thisFile, __LINE__, "" );
			muFail = mu;
			muNext = (muFail+muSafe)/2.0;
		end
	end
	%
	%
	% Do one tight iteration.
	mu = muCrit*( 1.0 + 2.0*opts.tol );
	matM = matH + (mu*matI);
	[ matR, cholFlag ] = chol( matM );
	if ( 0 ~= cholFlag )
		error( "chol() failed for penultimate." );
	end
	clear opts;
	opts.v0 = vecPsi; % Finally use vecPsi instead of vecG.
	opts.cholB = true;
	[ vecPsi, x, eigsFlag ] = eigs( matI, matR, 1, 'lm', opts );
	if ( 0 ~= eigsFlag )
		error( "eigs() failed for penultimate mu." );
	end
	muCrit = mu - (1.0./x);
	%
	lambdaSA = -muCrit;
return;
end
