clear
mydefs;
setprngstates( 0 );
%
sizeX = 5
numPts = 11
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
%%%matG += 0.1*randn(size(matG));

[ f_minWas, pt0 ] = min(rvecF);
%
prm = [];
prm.pt0 = pt0;
msg( __FILE__, __LINE__, "Calling hessfitss()..." );
tic();
[ matV, fss, vecGss, matHss, hessfitssDat ] = hessfitss( sizeX, numPts, matX, rvecF, matG, pt0, prm );
eigsOfMatHss = eig(matHss)
toc();
vecYNewton = -(matHss\vecGss);
vecXNext = matX(:,pt0) + matV * vecYNewton;
secret_vecDX0 = vecXNext - secret_vecX0;
f_minWas = f_minWas
vecYTest = -0.001*vecGss;
f_test = fss + ( vecGss'*vecYTest ) + ( vecYTest'*matHss*vecYTest )/2.0
f_expect = fss + ( vecGss'*vecYNewton ) + ( vecYNewton'*matHss*vecYNewton )/2.0
%vecNablaYF_expect = vecGss + matHss*vecYNewton
f_wouldBe = secret_f0 + (secret_vecDX0'*secret_matH*secret_vecDX0)/2.0
%
msg( __FILE__, __LINE__, "End of test." ); return;
