clear;
setprngstates(0);
numFigs = 0;
sizeX = 100
sizeF = 1
sizeK = ceil(sizeX^0.8)
sizeL = ceil(sizeK/2.0)
matV = utorthdrop( randn(sizeX,sizeK) );
%
matJSecret = zeros(sizeF,sizeX);
for nf=1:sizeF
for m=1:sizeL
	nx = 1 + floor(sizeX*(1.0-100.0*eps)*rand);
	nzeListSecret(m) = nx;
	matJSecret(nf,nx) += randn();
endfor
endfor
matW = matJSecret * matV;
%
if (0)
	prm = [];
	time0 = time(); [ matJApprox350, datOut ] = sja_scratch350( matV, matW, prm ); time350 = time()-time0;
	rd350 = reldiff( matJSecret, matJApprox350 );
	msg( __FILE__, __LINE__, sprintf( " sja_scratch350() produced %0.3e in %0.3f seconds.", rd350, time350 ) );
	return;
endif
%
prm = [];
time0 = time(); [ matJApprox300, datOut ] = sja_scratch300( matV, matW, prm ); time300 = time()-time0;
time0 = time(); [ matJApprox350, datOut ] = sja_scratch350( matV, matW, prm ); time350 = time()-time0;
rd300 = reldiff( matJSecret, matJApprox300 );
rd350 = reldiff( matJSecret, matJApprox350 );
msg( __FILE__, __LINE__, sprintf( " sja_scratch300() produced %0.3e in %0.3f seconds.", rd300, time300 ) );
msg( __FILE__, __LINE__, sprintf( " sja_scratch350() produced %0.3e in %0.3f seconds.", rd350, time350 ) );
return;
%
prm = [];
prm.maxNumNZEPerRow = floor( sizeK - sqrt(sizeK) );
time0 = time(); [ matJApprox100, datOut ] = sja_scratch100( matV, matW, prm ); time100 = time()-time0;
time0 = time(); [ matJApprox200, datOut ] = sja_scratch200( matV, matW, prm ); time200 = time()-time0;
time0 = time(); [ matJApprox300, datOut ] = sja_scratch300( matV, matW, prm ); time300 = time()-time0;
sort(nzeListSecret)
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
