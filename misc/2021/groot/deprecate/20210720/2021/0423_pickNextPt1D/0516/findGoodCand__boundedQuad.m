function [ xCand, meritCand ] = findGoodCand__boundedQuad( ...
  xl, xc, xr, gl, gc, gr, prm = [] );
  	% Should-be-precompiled...
	commondefs;
	thisFile = "findGoodCand__boundedQuad.m";
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	if (0)
		echo__thisFile = thisFile
		echo__xl = xl
		echo__xc = xc
		echo__xr = xr
		echo__gl = gl
		echo__gc = gc
		echo__gr = gr
	end
	msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
	  "{ %f, %f, %f }, { %f, %f, %f }.",xl,xc,xr,gl,gc,gr ) );
	%
	assert( isrealscalar(xl) );
	assert( isrealscalar(xc) );
	assert( isrealscalar(xr) );
	assert( isrealscalar(gl) );
	assert( isrealscalar(gc) );
	assert( isrealscalar(gr) );
	assert( xl < xc );
	assert( xc < xr );
	assert( gl > gc ); % Require monotonicity in g?
	assert( gc > gr );
	assert( gl > 0.0 );
	assert( gc != 0.0 );
	assert( gr < 0.0 );
	%
	%
	% Find merit.
	deltaXRatio = (xc-xl)/(xr-xc);
	deltaGRatio = (gc-gl)/(gr-gc);
	assert( deltaXRatio > 0.0 );
	assert( deltaGRatio > 0.0 );
	ldelta = log(deltaXRatio);
	lslope = log(deltaGRatio/deltaXRatio);
	meritExpCoeff = mygetfield( prm, "meritExpCoeff", 0.1 );
	assert( isrealscalar(meritExpCoeff) );
	assert( 0.0 < meritExpCoeff );
	meritCand = exp( -meritExpCoeff*( ldelta^2 + lslope^2 ) );
	%
	%
	% Find model root.
	% Normalize in x.
	yl = (xl-xc)/(xr-xl);
	yc = 0.0;
	yr = (xr-xc)/(xr-xl);
	yVec = [ yl; yc; yr ];
	gVec = [ gl; gc; gr ];
	yMat = [ ones(3,1), yVec, yVec.^2 ];
	cVec = yMat \ gVec;
	assert( isrealarray(cVec,[3,1]) );
	c0 = cVec(1);
	c1 = cVec(2);
	c2 = cVec(3);
	%
	% Watch out for co-linear case.
	if ( abs(4.0*c0*c2) <= eps*(c1^2) )
		assert( abs(c1) > eps*abs(c0) );
		yCand = -c0 / c1;
		xCand = xc + yCand*(xr-xl);
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "merit = %f, xCand = %f",meritCand,xCand ) );
		return;
	end
	%
	% There should always be exactly one root in the appropriate interval,
	% and it should correspond to the minus sign.
	discr = c1^2 - 4.0*c0*c2;
	assert( 0.0 <= discr );
	yCand = ( -c1 - sqrt(discr) ) / ( 2.0*c2 );
	assert( isrealscalar(yCand) );
	xCand =	xc + yCand*(xr-xl);
	msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
	  "merit = %f, xCand = %f",meritCand,xCand ) );
return;
end
