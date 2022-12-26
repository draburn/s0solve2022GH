function [ vecX, datOut ] = newtish_eig( matA, vecB, prm=[] )
	vecX = [];
	datOut = [];
	%
	assert( 2 <= nargin );
	assert( nargin <= 3 );
	sz = size(vecB,1);
	assert( isrealarray(vecB,[sz,1]) );
	assert( isrealarray(matA,[sz,sz]) );
	debugMode = mygetfield( prm, "debugMode", false );
	assert( issymmetric( matA ) );
	%
	[ matPsi, matLambda ] = eig( matA );
	assert( isrealarray(matLambda,[sz,sz]) );
	vecLambda = diag(matLambda);
	epsLambdaMin = mygetfield( prm, "epsLambdaMin", 1.0e-4 );
	lambdaMinImposed = epsLambdaMin * max(abs(vecLambda));
	assert( 0.0 < lambdaMinImposed );
	vecLambdaMod = vecLambda;
	vecLambdaMod( vecLambdaMod < lambdaMinImposed ) = lambdaMinImposed;
	%
	vecX = matPsi * ( (matPsi'*vecB) ./ vecLambdaMod );
	%
	datOut.matPsi = matPsi;
	datOut.vecLambda = vecLambda;
	datOut.vecLambdaMod = vecLambdaMod;
return;
endfunction

%!test
%!	setprngstates(0);
%!	%sizeX = 2 + ceil(20*rand());
%!	sizeX = 5
%!	size1 = sizeX-1;
%!	size2 = 1;
%!	foo = randn(size1,sizeX); matA1 = foo'*foo; clear foo;
%!	foo = randn(size2,sizeX); matA2 = foo'*foo; clear foo;
%!	vecX = randn(sizeX,1);
%!	vecBMod = randn(sizeX,1);
%!	epsA = 1.0e-6;
%!	%
%!	msg( __FILE__, __LINE__, "Strongly positive-definite test..." );
%!	matA = matA1 + 1.0*matA2;
%!	vecB = matA * vecX;
%!	vecXCalc = newtish_eig( matA, vecB );
%!	vecBCalc = matA * vecXCalc;
%!	assert( reldiff(vecB,vecBCalc) < eps^0.3 );
%!	assert( reldiff(vecX,vecXCalc) < eps^0.3 );
%!	%
%!	msg( __FILE__, __LINE__, "Weakly positive-definite test..." );
%!	matA = matA1 + epsA*matA2;
%!	vecB = matA * vecX;
%!	vecXCalc = newtish_eig( matA, vecB );
%!	assert( ~isempty(vecXCalc) );
%!	vecBCalc = matA * vecX;
%!	assert( reldiff(vecB,vecBCalc) < eps^0.3 );
%!	%
%!	msg( __FILE__, __LINE__, "Convergent positive-semi-definite test..." );
%!	matA = matA1 + 0.0*matA2;
%!	vecB = matA * vecX;
%!	vecXCalc = newtish_eig( matA, vecB );
%!	vecBCalc = matA * vecXCalc;
%!	assert( ~isempty(vecXCalc) );
%!	assert( reldiff(vecB,vecBCalc) < eps^0.3 );
%!	%
%!	msg( __FILE__, __LINE__, "Non-convergent positive-semi-definite test..." );
%!	matA = matA1 + 0.0*matA2;
%!	vecB = matA * vecX + 0.1*vecBMod;
%!	vecXCalc = newtish_eig( matA, vecB );
%!	vecBCalc = matA * vecXCalc;
%!	assert( ~isempty(vecXCalc) );
%!	%
%!	msg( __FILE__, __LINE__, "Strongly 'has negative' test..." );
%!	matA = -( matA1 + matA2 );
%!	vecB = matA * vecX;
%!	vecXCalc = newtish_eig( matA, vecB, false );
%!	vecBCalc = matA * vecXCalc;
%!	assert( ~isempty(vecXCalc) );
