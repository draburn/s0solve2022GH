clear;
mydefs;
setprngstates(15160912);
numFigs = 0;
sizeX = 100
sizeF = 1
sizeK = ceil(sizeX^0.8)
sizeL = ceil(sizeK/2.0)
matV = randn(sizeX,sizeK);
matV = utorthdrop( matV );
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
prm = [];
prm.maxNumNZEPerRow = floor( sizeK - sqrt(sizeK) );
%prm.verbLev = VERBLEV__INFO;
%sorted_nzeListSecret = sort(nzeListSecret)
%return;
%
if (1)
	time0 = time(); [ matJA_basic, datOut ] = sja_basic( matV, matW, prm ); time_basic= time()-time0;
	rd_basic = reldiff( matJSecret, matJA_basic );
	msg( __FILE__, __LINE__, sprintf( " sja_scratch100() produced %0.3e in %0.3f seconds.", rd_basic, time_basic ) );
	return;
endif
%
if (1)
	time0 = time(); [ matJApprox050, datOut ] = sja_scratch050( matV, matW, prm ); time050 = time()-time0;
	rd050 = reldiff( matJSecret, matJApprox050 );
	msg( __FILE__, __LINE__, sprintf( " sja_scratch050() produced %0.3e in %0.3f seconds.", rd050, time050 ) );
endif
%
if (0)
	time0 = time(); [ matJApprox100, datOut ] = sja_scratch100( matV, matW, prm ); time100 = time()-time0;
	rd100 = reldiff( matJSecret, matJApprox100 );
	msg( __FILE__, __LINE__, sprintf( " sja_scratch100() produced %0.3e in %0.3f seconds.", rd100, time100 ) );
endif
%
if (0)
	time0 = time(); [ matJApprox200, datOut ] = sja_scratch200( matV, matW, prm ); time200 = time()-time0;
	rd200 = reldiff( matJSecret, matJApprox200 );
	msg( __FILE__, __LINE__, sprintf( " sja_scratch100() produced %0.3e in %0.3f seconds.", rd100, time100 ) );
endif
%
if (0)
	time0 = time(); [ matJApprox300, datOut ] = sja_scratch300( matV, matW, prm ); time300 = time()-time0;
	rd300 = reldiff( matJSecret, matJApprox300 );
	msg( __FILE__, __LINE__, sprintf( " sja_scratch300() produced %0.3e in %0.3f seconds.", rd300, time300 ) );
endif
%
if (1)
	time0 = time(); [ matJApprox350, datOut ] = sja_scratch350( matV, matW, prm ); time350 = time()-time0;
	rd350 = reldiff( matJSecret, matJApprox350 );
	msg( __FILE__, __LINE__, sprintf( " sja_scratch350() produced %0.3e in %0.3f seconds.", rd350, time350 ) );
endif
%
%sorted_nzeListSecret = sort(nzeListSecret)
msg( __FILE__, __LINE__, "Goodbye." ); return;
%
assert( reldiff( matJApprox300, matJApprox350 ) < sqrt(eps) );
%
numFigs++; figure(numFigs);
plot( ...
  sumsq((matJApprox050-matJSecret)'), 'v-', ...
  sumsq((matJApprox100-matJSecret)'), 'x-', ...
  sumsq((matJApprox200-matJSecret)'), 's-', ...
  sumsq((matJApprox300-matJSecret)'), '^-' );
grid on;
%
numFigs++; figure(numFigs);
n = 1;
plot( ...
  matJApprox100(n,:), 'v-', ...
  matJApprox100(n,:), 'x-', ...
  matJApprox200(n,:), 's-', ...
  matJApprox300(n,:), '^-', ...
  matJSecret(n,:), 'o-' );
grid on;
