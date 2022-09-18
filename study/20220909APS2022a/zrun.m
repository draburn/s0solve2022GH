clear
mydefs;
printf("\n\n");
%
% SET PRM.
seedTime = mod( round(now*1E11), 1E8 ); % In case you want this.
probSetPrm = [];
probSetPrm.probType = "test1";
probSetPrm.numProbs = 10;
probSetPrm.numUnknowns = 10;
probSetPrm.setSeed = 0;
%
n = 0;
n++; a.s(n).f = @groot_fsolve;
n++; a.s(n).f = @groot_jfnk_basic;
n++; a.s(n).f = @groot_jfnk_convent;
n++; a.s(n).f = @groot_jfnk_baseline;
n++; a.s(n).f = @groot_jfnk_scratch100;
%
% Note: zapsrun__start converts "a" to "algoSetPrm".
zrun__start;
%
% DO WORK.
zcdo = zcomp( probSetPrm, algoSetPrm, default_solverPrm, zcompPrm );
%
zrun__finish; % Does nothing?
%
msg( __FILE__, __LINE__, "Goodbye." );
printf("\n\n");
