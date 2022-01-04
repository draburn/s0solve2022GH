clear;
commondefs;
thisFile = "test_fminunc";
setprngstates(0);
numFigs = 0;
%
sizeX = 10;
sizeF = sizeX;
vecXRoot = randn(sizeX,1);
matA = randn(sizeF,sizeX);
funchBigL = @(x)( 0.5*(norm(matA*(x-vecXRoot))^2) );
funchVecG = @(x)( matA'*matA*(x-vecXRoot) );
vecX0 = randn(sizeX,1);
%
fcn = @(x)({ funchBigL(x), funchVecG(x) });
fcn0 = fcn(vecX0)
opt = [];
opt.GradObj = "on";
%
msg( thisFile, __LINE__, "Calling findLocMinLossG_fminunc()..." );
%%%[ vecXF, retCode, datOut ] = findLocMinLossG_fminunc( funchBigL, funchVecG, vecX0 );
%%%fcn = @(x)([ funchBigL(x), funchVecG(x)' ]);
%%%fcn = @(x)([ funchBigL(x)  ]);
time0 = time();
vecXF = fminunc( fcn, vecX0 )
retCode = RETCODE__SUCCESS;
timeF = time();
%
bigL0 = funchBigL(vecX0);
bigLF = funchBigL(vecXF);
magD0 = norm(vecX0-vecXRoot);
magDF = norm(vecXF-vecXRoot);
msg( thisFile, __LINE__, sprintf( "findLocMinLossG_fminunc() returned retCode %s.", retcode2str(retCode) ) );
%msg( thisFile, __LINE__, sprintf( "  Performed %d iterations", datOut.iterCount ) );
msg( thisFile, __LINE__, sprintf( "  Loss went from %0.3e to %0.3e.", bigL0, bigLF ) );
msg( thisFile, __LINE__, sprintf( "  ||x-xRoot|| went from %0.3e to %0.3e.",magD0, magDF ) );
msg( thisFile, __LINE__, sprintf( "  Elapsed time was %gs.", time()-time0 ) );
