clear;
setprngstates(0);
%sizeK = 3; sizeX = 5;
sizeK = 4; sizeX = 100;
sizeF = sizeX;
%
matJSecret = diag( (1:sizeX) );
for nf=1:sizeF
	nx = 1 + floor(sizeX*(1.0-100.0*eps)*rand);
	matJSecret(nf,nx) += randn();
endfor
matV = utorthdrop( randn(sizeX,sizeK) );
matW = matJSecret * matV;
%
prm.maxNZEPerRow = 2;
[ matJApprox100, datOut ] = sja_scratch100( matV, matW, prm );
rd100 = reldiff( matJSecret, matJApprox100 )
[ matJApprox200, datOut ] = sja_scratch200( matV, matW, prm );
rd200 = reldiff( matJSecret, matJApprox200 )
