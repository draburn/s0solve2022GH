clear
mydefs;
setprngstates( 0 );
%
sizeX = 10
numPts = 51
%numPts = 200
numUnk = (sizeX*(sizeX+1))/2 + sizeX + 1
numGvn = (1+sizeX)*numPts
%
secret_f0 = 0.0;
secret_vecX0 = randn(sizeX,1);
foo = randn(sizeX,sizeX);
secret_matH = foo'*foo;
clear foo;
%
matX = randn(sizeX,numPts);
matX(:,1) = 0.0;
secret_matDX = matX - secret_vecX0;
rvecF = secret_f0 + sum( secret_matDX .* ( secret_matH * secret_matDX ), 1 )/2.0;
matG = secret_matH * secret_matDX;
%
matG += 0.1*randn(size(matG));
%
prm = [];
msg( __FILE__, __LINE__, "Calling hessfitss()..." );
tic();
[ matV, fss, vecGss, matHss, hessfitssDat ] = hessfitss( sizeX, numPts, matX, rvecF, matG, prm );
toc();
calc_vecX0ss = matHss\vecGss;
calc_vecX0 = -matV*calc_vecX0ss;
secret_vecDX0 = calc_vecX0 - secret_vecX0;
f_minWas = min(rvecF)
f_expect = fss + (vecGss'*calc_vecX0ss) + (calc_vecX0ss'*matHss*calc_vecX0ss)/2.0
f_wouldBe = secret_f0 + (secret_vecDX0'*secret_matH*secret_vecDX0)/2.0
%
msg( __FILE__, __LINE__, "End of test." ); return;
