clear;
setprngstates(0);
%sizeK = 3; sizeX = 5;
sizeK = 10; sizeX = 100;
sizeF = sizeX;
%
%%%matJSecret = diag( (1:sizeX) );
matJSecret = eye(sizeX);
for nf=1:sizeF
	nx = 1 + floor(sizeX*(1.0-100.0*eps)*rand);
	matJSecret(nf,nx) += randn();
endfor
matV = utorthdrop( randn(sizeX,sizeK) );
matW = matJSecret * matV;
%
prm = [];
prm.maxNumNZEPerRow = 9;
tic();
[ matJApprox100, datOut ] = sja_scratch100( matV, matW, prm );
toc();
tic();
[ matJApprox200, datOut ] = sja_scratch200( matV, matW, prm );
toc();
rd100 = reldiff( matJSecret, matJApprox100 )
rd200 = reldiff( matJSecret, matJApprox200 )
