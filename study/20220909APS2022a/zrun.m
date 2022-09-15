clear
mydefs;
printf("\n\n");
%
% SET PRM.
seedTime = mod( round(now*1E11), 1E8 ); % In case you want this.
probSetPrm = [];
probSetPrm.numProbs = 10;
probSetPrm.probType = "test1";
probSetPrm.numUnknowns = 20;
probSetPrm.setSeed = 0;
%
n = 0;
n++; a.s(n).f = @groot_fsolve;
n++; a.s(n).f = @groot_jfnk_basic;
n++; a.s(n).f = @groot_jfnk_basic; a.s(n).p.btCoeff = 0.0;
%
% Note: zapsrun__start converts "a" to "algoSetPrm".
zrun__start;
%
% DO WORK.
csdo = compsolv( probSetPrm, algoSetPrm, default_solverPrm, compsolvPrm );
%
%
zrun__finish;
msg( __FILE__, __LINE__, "Goodbye." );
printf("\n\n");
