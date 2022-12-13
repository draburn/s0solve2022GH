clear;
mydefs;
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
setprngstates(0);
sizeX = 5;
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
msg( __FILE__, __LINE__, "Performing positive-definite, test..." );
matA = matA1 + epsA*matA2;
vecB = matA * vecX;
vecXCalc = mycholdiv( matA, vecB, true );
vecBCalc = matA * vecXCalc;
xResults = [ norm(vecX), norm(vecXCalc), reldiff(vecX,vecXCalc) ]
bResults = [ norm(vecB), norm(vecBCalc), reldiff(vecB,vecBCalc) ]
deltaF = (vecXCalc'*matA*vecXCalc)/2.0 - vecXCalc'*vecB
vecXCalcPD = vecXCalc;
deltaFPD = deltaF;
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
msg( __FILE__, __LINE__, "Performing 'has negative', test..." );
matA = matA1 - epsA*matA2;
vecB = matA * vecX;
vecXCalc = mycholdiv( matA, vecB, false );
vecBCalc = matA * vecXCalc;
xResults = [ norm(vecX), norm(vecXCalc), reldiff(vecX,vecXCalc) ]
bResults = [ norm(vecB), norm(vecBCalc), reldiff(vecB,vecBCalc) ]
deltaF = (vecXCalc'*matA*vecXCalc)/2.0 - vecXCalc'*vecB
vecXCalcHN = vecXCalc;
deltaFHN = deltaF;
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
msg( __FILE__, __LINE__, "" );
%
matXCompare = [ vecXCalcPD, vecXCalcPSD, vecXCalcHN ]
rvecFCompare = [ deltaFPD, deltaFPSD, deltaFHN ]
%
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "" );
