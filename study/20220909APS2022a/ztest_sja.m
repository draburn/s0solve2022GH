clear;
setprngstates(0);
sizeX = 10;
sizeF = sizeX;
sizeK = 3;
%
matJ_true = diag( (1:sizeX) )
matV = utorthdrop( randn(sizeX,sizeK) )
matW = matJ_true * matV
%
prm.maxNumLEPerRow = 1;
[ matJApprox, datOut ] = sja_scratch000( matV, matW, prm )
