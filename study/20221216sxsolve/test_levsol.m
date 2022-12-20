clear;
setprngstates(0);
sizeX = 5;

msg( __FILE__, __LINE__, "~ Pos-def ~" );
fCrit = 0.0;
matH = mtm(randn(sizeX,sizeX));
vecXCrit = randn(sizeX,1);
funchD = @(x)( x - vecXCrit );
funchGOfD = @(d)( matH * d );
funchFOfD = @(d)( fCrit + sum( d .* (matH*d), 1 )/2.0 );
funchG = @(x)(funchGOfD(funchD(x)));
funchF = @(x)(funchFOfD(funchD(x)));
vecX0 = zeros(sizeX,1);
vecG0 = funchG(vecX0);
f0 = funchF(vecX0)

matB = diag((1:sizeX));

msg( __FILE__, __LINE__, "No scaling nor bMax..." );
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
%
msg( __FILE__, __LINE__, "No scaling nor bMax, fMin = 1.0..." );
prm = [];
prm.fMin = 1.0;
[ vecDelta, datOut ] = levsol_eig( f0, vecG0, matH, [], [], prm );
vecX1 = vecX0 + vecDelta;
b = norm(matB*vecDelta)
f1 = funchF(vecX1)


msg( __FILE__, __LINE__, "~ Likely has-neg ~" );
fCrit = 100.0;
matH = mtm(randn(sizeX,sizeX)) - mtm(randn(sizeX,sizeX));
minEig = min(eig(matH))
vecXCrit = randn(sizeX,1);
funchD = @(x)( x - vecXCrit );
funchGOfD = @(d)( matH * d );
funchFOfD = @(d)( fCrit + sum( d .* (matH*d), 1 )/2.0 );
funchG = @(x)(funchGOfD(funchD(x)));
funchF = @(x)(funchFOfD(funchD(x)));
vecX0 = zeros(sizeX,1);
vecG0 = funchG(vecX0);
f0 = funchF(vecX0)

msg( __FILE__, __LINE__, "No scaling nor bMax..." );
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



msg( __FILE__, __LINE__, "~ Likely pos-semi-def ~" );
fCrit = 0.0;
matH = mtm(randn(sizeX-1,sizeX));
minEig = min(eig(matH))
vecXCrit = randn(sizeX,1);
funchD = @(x)( x - vecXCrit );
funchGOfD = @(d)( matH * d );
funchFOfD = @(d)( fCrit + sum( d .* (matH*d), 1 )/2.0 );
funchG = @(x)(funchGOfD(funchD(x)));
funchF = @(x)(funchFOfD(funchD(x)));
vecX0 = zeros(sizeX,1);
vecG0 = funchG(vecX0);
f0 = funchF(vecX0)

msg( __FILE__, __LINE__, "No scaling nor bMax..." );
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
