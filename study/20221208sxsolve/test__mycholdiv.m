clear;
mydefs;
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
setprngstates(0);
sizeX = 10;
size1 = sizeX-1;
size2 = 1;
foo = randn(size1,sizeX); matA1 = foo'*foo; clear foo;
foo = randn(size2,sizeX); matA2 = foo'*foo; clear foo;
vecX = randn(sizeX,1);
vecBMod = randn(sizeX,1);
epsA = 1.0e-6;
%
%
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
msg( __FILE__, __LINE__, "Performing strongly positive-definite, test..." );
matA = matA1 + 1.0*matA2;
vecB = matA * vecX;
vecXCalc = mycholdiv( matA, vecB, true );
vecBCalc = matA * vecXCalc;
xResults = [ norm(vecX), norm(vecXCalc), reldiff(vecX,vecXCalc) ]
bResults = [ norm(vecB), norm(vecBCalc), reldiff(vecB,vecBCalc) ]
deltaF = (vecXCalc'*matA*vecXCalc)/2.0 - vecXCalc'*vecB
assert( reldiff(vecX,vecXCalc) < eps^0.3 );
assert( reldiff(vecB,vecBCalc) < eps^0.3 );
vecXCalcStrongPD = vecXCalc;
deltaFStrongPD = deltaF;
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
msg( __FILE__, __LINE__, "" );
%
%
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
msg( __FILE__, __LINE__, "Performing marginally positive-definite, test..." );
matA = matA1 + epsA*matA2;
vecB = matA * vecX;
vecXCalc = mycholdiv( matA, vecB, true );
vecBCalc = matA * vecXCalc;
xResults = [ norm(vecX), norm(vecXCalc), reldiff(vecX,vecXCalc) ]
bResults = [ norm(vecB), norm(vecBCalc), reldiff(vecB,vecBCalc) ]
deltaF = (vecXCalc'*matA*vecXCalc)/2.0 - vecXCalc'*vecB
assert( reldiff(vecX,vecXCalc) < eps^0.3 );
assert( reldiff(vecB,vecBCalc) < eps^0.3 );
vecXCalcMarginalPD = vecXCalc;
deltaFMarginalPD = deltaF;
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
msg( __FILE__, __LINE__, "" );
%
%
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
msg( __FILE__, __LINE__, "Performing convergant positive-semi-definite, test..." );
matA = matA1 + 0.0*matA2;
vecB = matA * vecX;
vecXCalc = mycholdiv( matA, vecB, true );
vecBCalc = matA * vecXCalc;
xResults = [ norm(vecX), norm(vecXCalc), reldiff(vecX,vecXCalc) ]
bResults = [ norm(vecB), norm(vecBCalc), reldiff(vecB,vecBCalc) ]
deltaF = (vecXCalc'*matA*vecXCalc)/2.0 - vecXCalc'*vecB
assert( reldiff(vecB,vecBCalc) < eps^0.3 );
vecXCalcPSD = vecXCalc;
deltaFPSD = deltaF;
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
msg( __FILE__, __LINE__, "" );
%
%
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
msg( __FILE__, __LINE__, "Performing non-convergant positive-semi-definite without-validation, test..." );
matA = matA1 + 0.0*matA2;
vecB = matA * vecX + 0.1*vecBMod;
vecXCalc = mycholdiv( matA, vecB, false );
vecBCalc = matA * vecXCalc;
xResults = [ norm(vecX), norm(vecXCalc), reldiff(vecX,vecXCalc) ]
bResults = [ norm(vecB), norm(vecBCalc), reldiff(vecB,vecBCalc) ]
deltaF = (vecXCalc'*matA*vecXCalc)/2.0 - vecXCalc'*vecB
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
msg( __FILE__, __LINE__, "" );
%
%
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
msg( __FILE__, __LINE__, "Performing non-convergant positive-semi-definite without-validation, test..." );
prm = [];
prm.validateExtrapolation = false;
vecXCalc = mycholdiv( matA, vecB, false, prm );
vecBCalc = matA * vecXCalc;
xResults = [ norm(vecX), norm(vecXCalc), reldiff(vecX,vecXCalc) ]
bResults = [ norm(vecB), norm(vecBCalc), reldiff(vecB,vecBCalc) ]
deltaF = (vecXCalc'*matA*vecXCalc)/2.0 - vecXCalc'*vecB
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
msg( __FILE__, __LINE__, "" );
%
%
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
msg( __FILE__, __LINE__, "Performing marginal 'has negative', test..." );
matA = matA1 - epsA*matA2;
vecB = matA * vecX;
vecXCalc = mycholdiv( matA, vecB, false );
vecBCalc = matA * vecXCalc;
xResults = [ norm(vecX), norm(vecXCalc), reldiff(vecX,vecXCalc) ]
bResults = [ norm(vecB), norm(vecBCalc), reldiff(vecB,vecBCalc) ]
deltaF = (vecXCalc'*matA*vecXCalc)/2.0 - vecXCalc'*vecB
vecXCalcMarginalHN = vecXCalc;
deltaFMarginalHN = deltaF;
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
msg( __FILE__, __LINE__, "" );
%
%
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
msg( __FILE__, __LINE__, "Performing strongly 'has negative', test..." );
matA = -( matA1 + matA2 );
vecB = matA * vecX;
vecXCalc = mycholdiv( matA, vecB, false );
vecBCalc = matA * vecXCalc;
xResults = [ norm(vecX), norm(vecXCalc), reldiff(vecX,vecXCalc) ]
bResults = [ norm(vecB), norm(vecBCalc), reldiff(vecB,vecBCalc) ]
deltaF = (vecXCalc'*matA*vecXCalc)/2.0 - vecXCalc'*vecB
vecXCalcStrongHN = vecXCalc;
deltaFStrongHN = deltaF;
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
msg( __FILE__, __LINE__, "" );
%
matXCompare = [ vecX, vecXCalcStrongPD, vecXCalcMarginalPD, vecXCalcPSD, vecXCalcMarginalHN, vecXCalcStrongHN ]
rvecFCompare = [ deltaFStrongPD, deltaFMarginalPD, deltaFPSD, deltaFMarginalHN, deltaFStrongHN ]
%
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "" );
