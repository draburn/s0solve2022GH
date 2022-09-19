clear;
setprngstates(0);
%sizeK = 3; sizeX = 5; sizeF = 5;
sizeK = 13; sizeX = 50; sizeF = 1;
matV = utorthdrop( randn(sizeX,sizeK) );
%
matJSecret = eye(sizeX);
for nf=1:sizeF
for m=1:3
	nx = 1 + floor(sizeX*(1.0-100.0*eps)*rand);
	matJSecret(nf,nx) += randn();
endfor
endfor
matJSecret = matJSecret(1,:);
matW = matJSecret * matV;
%
prm = [];
prm.maxNumNZEPerRow = sizeK-1;
tic();
[ matJApprox100, datOut ] = sja_scratch100( matV, matW, prm );
toc();
tic();
[ matJApprox200, datOut ] = sja_scratch200( matV, matW, prm );
toc();
rd100 = reldiff( matJSecret, matJApprox100 )
rd200 = reldiff( matJSecret, matJApprox200 )
plot( ...
  matJSecret, 'o-', ...
  matJApprox100, 'x-', ...
  matJApprox200, 's-' );
grid on;
