clear;
setprngstates(0);
%for n = [ 200, 500, 1000, 2000 ]
for n = [ 200, 500, 1000 ]
%for n = 5
	msg( __FILE__, __LINE__, sprintf( "Timing for n = %d...", n ) );
	f = randn(n,1);
	%
	t = time();
	j = randn(n,n);
	%%%j = randn(round(n/2),n);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"j = randn(%d,%d)\".", 1000.0*(time()-t), n, n) );
	%%%j .*= exp( 3.0*randn(n,n) );
	%
	j2_temp = randn(n,n);
	t = time();
	j3_temp = j + j2_temp;
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"j3 = j1 + j2\".", 1000.0*(time()-t)) );
	clear j2_temp;
	clear j3_temp;
	%
	x_temp = randn(n,1);
	t = time();
	f_temp = j * x_temp;
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"f = j * x\".", 1000.0*(time()-t)) );
	clear x_temp;
	clear f_temp;
	%
	t = time();
	h = j' * j;
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"h = j' * j\".", 1000.0*(time()-t)) );
	%%%h += norm(diag(h))*(eps^0.8)*eye(n,n);
	%%%h = diag(diag(j' * j));
	%
	t = time();
	r = chol(h);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"r = chol(h)\".", 1000.0*(time()-t)) );
	%
	t = time();
	x = r \ ( r' \ f );
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"x = r \\ ( r' \\ f )\".", 1000.0*(time()-t)) );
	%
	t = time();
	x = h \ f;
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"x = h \\ f\".", 1000.0*(time()-t)) );
	%
	%
	%
	t = time();
	[ l, u ] = lu(h);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"[ l, u ] = lu(h)\".", 1000.0*(time()-t)) )
	%
	t = time();
	hNorm = norm(h);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"hNorm = norm(h)\".", 1000.0*(time()-t)) );
	%
	t = time();
	lambda = eig(h);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"lambda = eig(h)\".", 1000.0*(time()-t)) );
	%
	t = time();
	[ matPsi, matLambda ] = eig(h);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"[ matPsi, matLambda ] = eig(h)\".", 1000.0*(time()-t)) );
	%
	t = time();
	[ matPsi, matLambda ] = eig(triu(h));
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"[ matPsi, matLambda ] = eig(triu(h))\".", 1000.0*(time()-t)) );
	%
	t = time();
	x = r \ ( r' \ h );
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"x = r \\ ( r' \\ h )\".", 1000.0*(time()-t)) );
	%
	t = time();
	hFro = norm(h,"fro");
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"hFro = norm(h,\"fro\")\".", 1000.0*(time()-t)) );
	%
	t = time();
	rc = rcond(h);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"rc = rcond(h)\".", 1000.0*(time()-t)) )
	%
	t = time();
	cEst = condest(h);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"cEst = condest(h)\".", 1000.0*(time()-t)) );
	%
	t = time();
	c = cond(h);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"c = cond(h)\".", 1000.0*(time()-t)) );
	%
	if ( n <= 500 )
		%
		t = time();
		[ v, lambda ] = eig(h);
		msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"[ v, lambda ] = eig(h)\".", 1000.0*(time()-t)) );
		%
		t = time();
		hInv = pinv(h);
		msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"hInv = pinv(h)\".", 1000.0*(time()-t)) );
		%
		t = time();
		q = orth(j);
		msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"q = orth(j)\".", 1000.0*(time()-t)) );
		%
		t = time();
		qtj = q' * j;
		msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"qtj = q' * j\".", 1000.0*(time()-t)) );
		%
		t = time();
		q = utorthdrop(j);
		msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"q = utorthdrop(j)\".", 1000.0*(time()-t)) );
		% Mine is faster!
		%
		t = time();
		phi = null(j);
		msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"phi = null(j)\".", 1000.0*(time()-t)) );
		%
		t = time();
		phi = null(h);
		msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"phi = null(h)\".", 1000.0*(time()-t)) );
	endif
	%
	for p = 0.1.^[ 0, 1, 2, 3, 4 ]
		j = randn(n,n);
		h = j'*j;
		h += p*diag(diag(h));
		c = condest(h);
		f = randn(n,1);
		g = j'*f;
		%
		t = time();
		x = cgs( h, g, 1e-2, n );
		msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"x = cgs( h, g, 1e-2, n )\" to %10.3e for condest %10.3e.", ...
		  1000.0*(time()-t), norm(g-h*x)/norm(g), c ) );
		%
		if ( c < 1e4 )
			% My (unoptimized!) linear sovler is not great, but comprable to built-in.
			t = time();
			funchH = @(z)( h*z );
			prm = [];
			prm.tol = 1e-2;
			x = linsolf( funchH, g, zeros(n,1), prm );
			msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"x = linsolf( h, g, ., 1e-2 )\" to %10.3e for condest %10.3e.", ...
			  1000.0*(time()-t), norm(g-h*x)/norm(g), c ) );
			%
			t = time();
			x = gmres( h, g );
			msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"x = gmres( h, g )\" to %10.3e for condest %10.3e.", ...
			  1000.0*(time()-t), norm(g-h*x)/norm(g), c ) );
		endif
		%
		if (0)
			% These are slower and have annoying output.
			t = time();
			x = bicg( h, g );
			msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"x = bicg( h, g )\" to %10.3e for condest %10.3e.", ...
			  1000.0*(time()-t), norm(g-h*x)/norm(g), c ) );
			%
			t = time();
			x = bicgstab( h, g );
			msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"x = bicgstab( h, g )\" to %10.3e for condest %10.3e.", ...
			  1000.0*(time()-t), norm(g-h*x)/norm(g), c ) );
		endif
	endfor
	%
	%
	continue;
	%
	%
	% BELOW IS HACKY STUFF....
	%
	t = time();
	hScl = norm(diag(h));
	sDiag = sqrt( diag(h)/hScl + eps );
	sInvDiag = 1.0./sDiag;
	s = diag(sDiag);
	sInv = diag(sInvDiag);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for generating sInv.", 1000.0*(time()-t)) );
	%
	t = time();
	hHat = ( sInv * h * sInv ) / hScl;
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"hHat = ( sInv * h * sInv ) / hScl\".", 1000.0*(time()-t)) );
	%
	rBefore = chol(h);
	rAfter = chol(hHat);
	[ max(diag(rBefore))/min(diag(rBefore)), max(diag(rAfter))/min(diag(rAfter)) ]
	%h
	%hHat
	%[ cond(h), cond(hHat) ]
	%continue
	%
	%
	%
	if (0)
	t = time();
	phi = randn(n,1);
	phi /= norm(phi);
	eta = h*phi;
	hScl = norm(eta);
	sDiag = sqrt( abs(eta)/hScl + eps );
	sInvDiag = 1.0./sDiag;
	s = diag(sDiag);
	sInv = diag(sInvDiag);
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for generating alt sInv.", 1000.0*(time()-t)) );
	%
	t = time();
	hHat = ( sInv * h * sInv ) / hScl;
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for \"hHat = ( sInv * h * sInv ) / hScl\".", 1000.0*(time()-t)) );
	%
	h
	hHat
	[ cond(h), cond(hHat) ]
	endif
	%
	%
	%
	j = randn(n-1,n);
	h = j'*j;
	t = time();
	hiHat = orth(h*eye(n,n-1));
	phi = randn(n,1);
	phi -= hiHat*(hiHat'*phi);
	phi -= hiHat*(hiHat'*phi);
	assert( norm(phi) > 0.0 );
	phi /= norm(phi);
	assert( norm(h*phi) < sqrt(eps) );
	msg( __FILE__, __LINE__, sprintf( "  %10.3f ms for finding null space.", 1000.0*(time()-t)) );
endfor
