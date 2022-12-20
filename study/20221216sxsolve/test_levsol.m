clear;
setprngstates(0);
sizeX = 5;

fRoot = 0.0;
matH = mtm(randn(sizeX,sizeX));
vecXRoot = randn(sizeX,1);
funchD = @(x)( x - vecXRoot );
funchGOfD = @(d)( matH * d );
funchFOfD = @(d)( fRoot + sum( d .* (matH*d), 1 )/2.0 );
funchG = @(x)(funchGOfD(funchD(x)));
funchF = @(x)(funchFOfD(funchD(x)));

vecX0 = zeros(sizeX,1);
vecG0 = funchG(vecX0);
f0 = funchF(vecX0)
matB = diag((1:sizeX));

msg( __FILE__, __LINE__, "No scaling or bMax..." );
[ vecDelta, datOut ] = levsol_eig( f0, vecG0, matH );
vecX1 = vecX0 + vecDelta;
b = norm(vecDelta)
f1 = funchF(vecX1)

msg( __FILE__, __LINE__, "Scaled but no bMax..." );
[ vecDelta, datOut ] = levsol_eig( f0, vecG0, matH, matB, [] );
vecX1 = vecX0 + vecDelta;
b = norm(matB*vecDelta)
f1 = funchF(vecX1)
%
msg( __FILE__, __LINE__, "Scaled with bMax = 1.0..." );
[ vecDelta, datOut ] = levsol_eig( f0, vecG0, matH, matB, 1.0 );
vecX1 = vecX0 + vecDelta;
b = norm(matB*vecDelta)
f1 = funchF(vecX1)
%
msg( __FILE__, __LINE__, "Scaled with no bMax, fMin = 1.0..." );
prm = [];
prm.fMin = 1.0;
[ vecDelta, datOut ] = levsol_eig( f0, vecG0, matH, matB, [], prm );
vecX1 = vecX0 + vecDelta;
b = norm(matB*vecDelta)
f1 = funchF(vecX1)
