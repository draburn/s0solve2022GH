clear
mydefs;
msg( __FILE__, __LINE__, "Performing dog leg test..." );
%
seedTime = mod( round(now*1E11), 1E8 ); % In case you want this.
strProbType = "lintest1";
bigN0 = 50;
probSeed = 0;
genFunchFTPrm = [];
[ funchFOfX, vecX0, genFunchFTDatOut ] = genFunchAPS2022_fromType( strProbType, bigN0, probSeed, genFunchFTPrm );
%
solverPrm = [];
solverPrm.testDogLeg = true;
groot_jfnk_convent( funchFOfX, vecX0, solverPrm );
%
msg( __FILE__, __LINE__, "Goodbye." );
