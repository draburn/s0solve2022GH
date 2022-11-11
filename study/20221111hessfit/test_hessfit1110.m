clear
mydefs;
setprngstates( 0 );
%
sizeX = 5;
numPts = 20;
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
msg( __FILE__, __LINE__, "Begin test main..." );
sizeFGL = 1 + sizeX + ((sizeX*(sizeX+1))/2)
prm = [];
%prm.rvecW1 = ones(1,numPts);
[ calc_f, calc_vecG, calc_matH, datOut ] = hessfit1110( matX, rvecF, matG, prm );
%
msg( __FILE__, __LINE__, "End of test." );
