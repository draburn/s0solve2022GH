clear;
mydefs;
setprngstates(20578384);
numFigs = 0;
sizeX = 100
sizeF = 100
sizeK = 90
sizeL = 3%ceil(sizeK/2.0)
%
matV = randn(sizeX,sizeK);
matJSecret = eye(sizeF,sizeX);
for nf=1:sizeF
for m=1:sizeL-1
	nx = 1 + floor(sizeX*(1.0-100.0*eps)*rand);
	nzeListSecret(m) = nx;
	matJSecret(nf,nx) += randn();
endfor
endfor
%
%matJSecret
matV = utorthdrop( matV );
matW = matJSecret*matV;
%
prm = [];
%prm.maxNumNZEPerRow = 5;
%tic(); [ matJA, datOut ] = sja_scratch350( matV, matW, prm ); toc(); rd = reldiff( matJSecret, matJA )
tic(); [ matJA, datOut ] = sja_basic( matV, matW, prm ); toc(); if ( ~isempty(matJA)) rd = reldiff( matJSecret, matJA ) endif
tic(); [ matJA, datOut ] = sja_fast( matV, matW, prm ); toc(); rd = reldiff( matJSecret, matJA )
tic(); [ matJA, datOut ] = sja_corr( matV, matW, prm ); toc(); rd = reldiff( matJSecret, matJA )
