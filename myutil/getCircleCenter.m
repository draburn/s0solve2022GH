function vecXC = getCircleCenter( vecX1, vecX2, vecX3 )
	vecA = vecX1 - vecX3;
	vecB = vecX2 - vecX3;
	aa = vecA'*vecA;
	bb = vecB'*vecB;
	ab = vecA'*vecB;
	if ( abs(aa*bb-ab*ab) < (eps^0.5)*(aa+bb+abs(ab)) )
		warning( "Points appear to be colinear.");
		vecXC = [];
		return;
	end
	matM = [ aa, ab; ab, bb ];
	vecV = 0.5*[ aa; bb ];
	vecC = matM\vecV;
	vecXC = vecX3 + vecC(1)*vecA + vecC(2)*vecB;
return

%!test
%!	clear;
%!	commondefs;
%!	thisFile = "test getCircleCenter 1";
%!	setprngstates();
%!	numFigs = 0;
%!	tic();
%!	%
%!	sizeX = 10;
%!	secret_vecXC = randn(sizeX,1);
%!	secret_bigR = abs(randn);
%!	secret_vecD1 = randn(sizeX,1);
%!	secret_vecD2 = randn(sizeX,1);
%!	secret_vecD3 = randn(sizeX,1);
%!	%
%!	vecX1 = secret_vecXC + secret_bigR*secret_vecD1/norm(secret_vecD1);
%!	vecX2 = secret_vecXC + secret_bigR*secret_vecD2/norm(secret_vecD2);
%!	vecX3 = secret_vecXC + secret_bigR*secret_vecD3/norm(secret_vecD3);
%!	%
%!	vecXC = getCircleCenter( vecX1, vecX2, vecX3 );
%!	%
%!	bigR = norm(vecX1 - vecXC);
%!	assert( abs(norm(vecX1-vecXC)-bigR) < (eps^0.75)*bigR )
%!	assert( abs(norm(vecX2-vecXC)-bigR) < (eps^0.75)*bigR )
%!	assert( abs(norm(vecX3-vecXC)-bigR) < (eps^0.75)*bigR )

%!test
%!	clear;
%!	commondefs;
%!	thisFile = "test getCircleCenter 2";
%!	setprngstates();
%!	numFigs = 0;
%!	tic();
%!	%
%!	sizeX = 2;
%!	secret_vecXC = randn(sizeX,1);
%!	secret_bigR = abs(randn);
%!	secret_vecD1 = randn(sizeX,1);
%!	secret_vecD2 = randn(sizeX,1);
%!	%
%!	vecX1 = secret_vecXC + secret_bigR*secret_vecD1/norm(secret_vecD1);
%!	vecX2 = secret_vecXC + secret_bigR*secret_vecD2/norm(secret_vecD2);
%!	secret_u = randn();
%!	vecX3 = secret_u*vecX1 + (1.0-secret_u)*vecX2;
%!	%
%!	vecXC = getCircleCenter( vecX1, vecX2, vecX3 );
%!	%
%!	assert( isempty(vecXC) )
%!	return
