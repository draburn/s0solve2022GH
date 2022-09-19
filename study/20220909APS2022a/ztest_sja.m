clear;
setprngstates(81560416);
numFigs = 0;
%sizeK = 3; sizeX = 5; sizeF = 5;
%sizeK = 10; sizeX = 100; sizeF = 1;
sizeK = 20; sizeX = 100; sizeF = 1;
%sizeK = 30; sizeX = 100; sizeF = sizeX;
matV = utorthdrop( randn(sizeX,sizeK) );
%
matJSecret = zeros(sizeF,sizeX);
for nf=1:sizeF
for m=1:10
	nx = 1 + floor(sizeX*(1.0-100.0*eps)*rand)
	matJSecret(nf,nx) += randn();
endfor
endfor
matW = matJSecret * matV;
%
prm = [];
prm.maxNumNZEPerRow = floor( sizeK - sqrt(sizeK) );
time0 = time(); [ matJApprox100, datOut ] = sja_scratch100( matV, matW, prm ); time100 = time()-time0;
time0 = time(); [ matJApprox200, datOut ] = sja_scratch200( matV, matW, prm ); time200 = time()-time0;
time0 = time(); [ matJApprox300, datOut ] = sja_scratch300( matV, matW, prm ); time300 = time()-time0;
rd100 = reldiff( matJSecret, matJApprox100 );
rd200 = reldiff( matJSecret, matJApprox200 );
rd300 = reldiff( matJSecret, matJApprox300 );
msg( __FILE__, __LINE__, sprintf( " sja_scratch100() produced %0.3e in %0.3f seconds.", rd100, time100 ) );
msg( __FILE__, __LINE__, sprintf( " sja_scratch200() produced %0.3e in %0.3f seconds.", rd200, time200 ) );
msg( __FILE__, __LINE__, sprintf( " sja_scratch300() produced %0.3e in %0.3f seconds.", rd300, time300 ) );
%
numFigs++; figure(numFigs);
plot( ...
  sumsq((matJApprox100-matJSecret)'), 'x-', ...
  sumsq((matJApprox200-matJSecret)'), 's-', ...
  sumsq((matJApprox300-matJSecret)'), '^-' );
grid on;
%
numFigs++; figure(numFigs);
n = 1;
n = 1;
plot( ...
  matJApprox100(n,:), 'x-', ...
  matJApprox200(n,:), 's-', ...
  matJApprox300(n,:), '^-', ...
  matJSecret(n,:), 'o-' );
grid on;
