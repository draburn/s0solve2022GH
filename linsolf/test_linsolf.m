if (0)
clear;
commondefs;
thisFile = "test_linsolf";
msg( thisFile, __LINE__, "Performing WIP test..." );
probSize = 800;
randnstate_before = randn("state");
randn("state",468);
matA = eye(probSize,probSize) + (0.04 * randn(probSize,probSize));
vecXSecret = randn(probSize,1);
vecB = matA * vecXSecret;
randn("state",randnstate_before);
funchMatAProd = @(vecV)( matA * vecV );
prm.verbLev = VERBLEV__COPIOUS;
[ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB, prm );
msg( thisFile, __LINE__, sprintf("retCode = %s.", retcode2str(retCode)) );
assert( RETCODE__NOT_SET == retCode );
msg( thisFile, __LINE__, "Finished WIP test.\n" );
end
%
%
clear;
commondefs;
thisFile = "test_linsolf";
msg( thisFile, __LINE__, "Performing vecB = 0 test..." );
probSize = 800;
randnstate_before = randn("state");
randn("state",468);
matA = eye(probSize,probSize) + (0.04 * randn(probSize,probSize));
vecXSecret = zeros(probSize,1);
vecB = matA * vecXSecret;
randn("state",randnstate_before);
funchMatAProd = @(vecV)( matA * vecV );
[ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB );
assert( 0.0 == sum(vecX.^2) );
assert( RETCODE__SUCCESS == retCode );
msg( thisFile, __LINE__, "Finished vecB = 0 test.\n" );
%
%
clear;
commondefs;
thisFile = "test_linsolf";
msg( thisFile, __LINE__, "Performing basic test..." );
probSize = 800;
randnstate_before = randn("state");
randn("state",468);
matA = eye(probSize,probSize) + (0.04 * randn(probSize,probSize));
vecXSecret = randn(probSize,1);
vecB = matA * vecXSecret;
randn("state",randnstate_before);
funchMatAProd = @(vecV)( matA * vecV );
prm.fracResTol = 1e-2;
prm.verbLev = VERBLEV__PROGRESS;
prm.reportInterval = 0.1;
[ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB, prm );
assert( RETCODE__SUCCESS == retCode );
vecRho = (matA*vecX) - vecB;
rhoNorm = sqrt(sum(vecRho.^2));
bNorm = sqrt(sum(vecB.^2));
res = rhoNorm / bNorm;
assert( res <= prm.fracResTol );
msg( thisFile, __LINE__, "Finished basic test.\n" );
%
%
clear;
commondefs;
thisFile = "test_linsolf";
msg( thisFile, __LINE__, "Performing numIterLimit test..." );
probSize = 1000;
randnstate_before = randn("state");
randn("state",661);
matA = eye(probSize,probSize) + (0.04 * randn(probSize,probSize));
vecXSecret = randn(probSize,1);
vecB = matA * vecXSecret;
randn("state",randnstate_before);
funchMatAProd = @(vecV)( matA * vecV );
prm.verbLev = VERBLEV__PROGRESS;
prm.fracResTol = 1e-2;
prm.numIterLimit = 10;
[ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB, prm );
assert( RETCODE__IMPOSED_STOP == retCode );
vecRho = (matA*vecX) - vecB;
rhoNorm = sqrt(sum(vecRho.^2));
bNorm = sqrt(sum(vecB.^2));
res = rhoNorm / bNorm;
assert( isfield(datOut,"numIter") );
assert(  (datOut.numIter == prm.numIterLimit) || (res <= prm.fracResTol) );
msg( thisFile, __LINE__, "Finished numIterLimit test.\n" );
%
%
clear;
commondefs;
thisFile = "test_linsolf";
msg( thisFile, __LINE__, "Performing projected vector is zero test..." );
probSize = 10;
randnstate_before = randn("state");
randn("state",122);
matA = zeros(probSize,probSize);
matA(2,1) = 1.0;
matA(3,2) = 1.0;
vecB = zeros(probSize,1);
vecB(1) = 1.0;
randn("state",randnstate_before);
funchMatAProd = @(vecV)( matA * vecV );
prm = [];
prm.verbLev = VERBLEV__MAIN;
[ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB, prm );
assert( RETCODE__ALGORITHM_BREAKDOWN == retCode );
assert( isfield(datOut,"numIter") );
assert( 2 == datOut.numIter );
msg( thisFile, __LINE__, "Finished projected vector is zero test.\n" );
%
%
clear;
commondefs;
thisFile = "test_linsolf";
msg( thisFile, __LINE__, "Performing non-invertible matrix test..." );
probSize = 100;
randnstate_before = randn("state");
randn("state",287);
a0Size = 20;
matA0 = eye(a0Size,probSize) + (0.1 * randn(a0Size,probSize));
matA = matA0.' * matA0;
vecB = randn(probSize,1);
randn("state",randnstate_before);
funchMatAProd = @(vecV)( matA * vecV );
prm = [];
prm.verbLev = VERBLEV__PROGRESS;
[ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB, prm );
assert( RETCODE__ALGORITHM_BREAKDOWN == retCode );
assert( isfield(datOut,"numIter") );
assert( a0Size+1 <= datOut.numIter ); % Would be ==, but FP issues exist.
assert( 1.1*(a0Size+1) >= datOut.numIter ); % A reasonable upper bound???
msg( thisFile, __LINE__, "Finished non-invertible matrix test.\n" );
%
%
if (0)
clear;
commondefs;
thisFile = "test_linsolf";
msg( thisFile, __LINE__, "Performing large-N intrinsic speed test..." );
probSize = 10000000;
randnstate_before = randn("state");
randn("state",167);
matA = diag(1.0+(0.5*randn(probSize,1)));
vecXSecret = randn(probSize,1);
vecB = matA * vecXSecret;
randn("state",randnstate_before);
funchMatAProd = @(vecV)( matA * vecV );
prm.fracResTol = 5.0e-2;
prm.verbLev = VERBLEV__PROGRESS;
time0 = time();
[ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB, prm );
msg( thisFile, __LINE__, sprintf("Took %f seconds.", time()-time0) );
assert( RETCODE__SUCCESS == retCode );
vecRho = (matA*vecX) - vecB;
rhoNorm = sqrt(sum(vecRho.^2));
bNorm = sqrt(sum(vecB.^2));
res = rhoNorm / bNorm;
assert( res <= prm.fracResTol );
msg( thisFile, __LINE__, "Finished large-N internal speed test.\n" );
%
%
clear;
commondefs;
thisFile = "test_linsolf";
msg( thisFile, __LINE__, "Performing large-k intrinsic speed test..." );
probSize = 1000;
randnstate_before = randn("state");
randn("state",167);
matA = diag(randn(probSize,1));
vecXSecret = randn(probSize,1);
vecB = matA * vecXSecret;
randn("state",randnstate_before);
funchMatAProd = @(vecV)( matA * vecV );
prm.fracResTol = 1.0e-6;
prm.verbLev = VERBLEV__PROGRESS;
time0 = time();
[ vecX, retCode, datOut ] = linsolf( funchMatAProd, vecB, prm );
msg( thisFile, __LINE__, sprintf("Took %f seconds.", time()-time0) );
assert( RETCODE__SUCCESS == retCode );
vecRho = (matA*vecX) - vecB;
rhoNorm = sqrt(sum(vecRho.^2));
bNorm = sqrt(sum(vecB.^2));
res = rhoNorm / bNorm;
assert( res <= prm.fracResTol );
msg( thisFile, __LINE__, "Finished large-k internal speed test.\n" );
end
