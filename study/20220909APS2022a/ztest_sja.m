clear;
setprngstates(0);
sizeX = 5;
sizeF = sizeX;
sizeK = 3;
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
[ matJApprox, datOut ] = sja_scratch200( matV, matW, prm );
reldiff( matJSecret, matJApprox )
