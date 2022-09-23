clear
mydefs;
printf("\n\n");
%
% SET PRM.
seedTime = mod( round(now*1E11), 1E8 ); % In case you want this.
probSetPrm = [];
probSetPrm.probType = "sja600";
probSetPrm.numProbs = 10;
probSetPrm.numUnknowns = 200;
probSetPrm.setSeed = 0;
%
n = 0;
n++; a.s(n).f = @groot_fsolve;
n++; a.s(n).f = @groot_jfnk_basic;
n++; a.s(n).f = @groot_jfnk_convent;
n++; a.s(n).f = @groot_jfnk_baseline;
n++; a.s(n).f = @groot_jfnk_sja;
%
% Note: zapsrun__start converts "a" to "algoSetPrm".
zrun__start;
default_solverPrm.linsolfPrm.tol = 0.01;
%
% DO WORK.
zcdo = zcomp( probSetPrm, algoSetPrm, default_solverPrm, zcompPrm );
%
zrun__finish; % Does nothing?
%
msg( __FILE__, __LINE__, "Goodbye." );
printf("\n\n");
