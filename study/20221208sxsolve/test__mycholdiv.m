clear;
mydefs;
setprngstates(0);
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
%
msg( __FILE__, __LINE__, "Performing basic, positive-definite, test..." );
sizeX = 5;
sizeF = sizeX;
foo = randn(sizeF,sizeX); matA = foo'*foo; clear foo;
vecX = randn(sizeX,1);
vecB = matA * vecX;
normX = norm(vecX)
vecXCalc = mycholdiv( matA, vecB );
vecBCalc = matA*vecXCalc;
normXCalc = norm(vecXCalc)
reldiffX = reldiff( vecXCalc, vecX )
normB = norm(vecB)
normBCalc = norm(vecBCalc)
reldiffB = reldiff( vecBCalc, vecB )
diffF = 0.5*(vecXCalc'*matA*vecXCalc) - vecXCalc'*vecB
assert( reldiffX < sqrt(eps) );
assert( reldiffB < sqrt(eps) );
msg( __FILE__, __LINE__, "Completed basic, positive-definite, test." );
clear matA;
clear vecX;
clear vecB;
clear vecXCalc;
clear vecBCalc;
msg( __FILE__, __LINE__, "--------------------------------------------------" );
%
msg( __FILE__, __LINE__, "Performing convergent semi-positive-definite, test..." );
sizeX = 5;
sizeF = 1;
foo = randn(sizeF,sizeX); matA = foo'*foo; clear foo;
vecX = randn(sizeX,1);
vecB = matA * vecX;
vecXCalc = mycholdiv( matA, vecB );
vecBCalc = matA*vecXCalc;
normX = norm(vecX)
normXCalc = norm(vecXCalc)
normB = norm(vecB)
normBCalc = norm(vecBCalc)
reldiffX = reldiff( vecXCalc, vecX )
reldiffB = reldiff( vecBCalc, vecB )
diffF = 0.5*(vecXCalc'*matA*vecXCalc) - vecXCalc'*vecB
assert( reldiffB < sqrt(eps) );
msg( __FILE__, __LINE__, "Completed convergent semi-positive-definite test." );
clear matA;
clear vecX;
clear vecB;
clear vecXCalc;
clear vecBCalc;
msg( __FILE__, __LINE__, "--------------------------------------------------" );
%
msg( __FILE__, __LINE__, "Performing non-convergent semi-positive-definite test (NO SOLUTION!)..." );
sizeX = 5;
sizeF = sizeX-1;
foo = randn(sizeF,sizeX); matA = foo'*foo; clear foo;
vecB = randn(sizeX,1);
vecXCalc = mycholdiv( matA, vecB, true );
vecBCalc = matA*vecXCalc;
normXCalc = norm(vecXCalc)
normB = norm(vecB)
normBCalc = norm(vecBCalc)
reldiffB = reldiff( vecBCalc, vecB )
diffF = 0.5*(vecXCalc'*matA*vecXCalc) - vecXCalc'*vecB
msg( __FILE__, __LINE__, "Completed non-convergent semi-positive-definite test (NO SOLUTION!)." );
clear matA;
clear vecB;
clear vecXCalc;
clear vecBCalc;
msg( __FILE__, __LINE__, "--------------------------------------------------" );
%
msg( __FILE__, __LINE__, "Performing 'has negative' test (NO SOLUTION!)..." );
sizeX = 5;
sizeFP = 3;
sizeFM = 1;
fooP = randn(sizeFP,sizeX);
fooM = randn(sizeFM,sizeX);
matA = fooP'*fooP - fooM'*fooM;
clear fooP;
clear fooM;
vecB = randn(sizeX,1);
vecXCalc = mycholdiv( matA, vecB, true );
vecBCalc = matA*vecXCalc;
normXCalc = norm(vecXCalc)
normB = norm(vecB)
normBCalc = norm(vecBCalc)
reldiffB = reldiff( vecBCalc, vecB )
diffF = 0.5*(vecXCalc'*matA*vecXCalc) - vecXCalc'*vecB
msg( __FILE__, __LINE__, "Completed 'has negative' test (NO SOLUTION!)." );
clear matA;
clear vecB;
clear vecXCalc;
clear vecBCalc;
%msg( __FILE__, __LINE__, "--------------------------------------------------" );
%
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
