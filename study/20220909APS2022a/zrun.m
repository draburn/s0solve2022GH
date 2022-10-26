clear
mydefs;
printf("\n\n");
%
% SET PRM.
seedTime = mod( round(now*1E11), 1E8 ); % In case you want this.
probSetPrm = [];
probSetPrm.probType = "aps2022base";
probSetPrm.numProbs = 20;
probSetPrm.numUnknowns = 100;
probSetPrm.setSeed = 0;
%
n = 0;
n++; a.s(n).f = @groot_fsolve;
%n++; a.s(n).f = @groot_jfnk_basic;
%n++; a.s(n).f = @groot_jfnk_convent;
n++; a.s(n).f = @groot_jfnk_baseline;
n++; a.s(n).f = @groot_jfnk_sja;
n++; a.s(n).f = @groot_jfnk_sja_looptr;
%n++; a.s(n).f = @groot_findZero_800;
%n++; a.s(n).f = @groot_jfnk_sja; a.s(n).p.linsolfPrm.sjaMethod = "omp";
%
% Note: zapsrun__start converts "a" to "algoSetPrm".
zrun__start;
default_solverPrm.linsolfPrm.tol = 0.01;
default_solverPrm.fevalLimit = 1E4;
%
% DO WORK.
%zcompPrm.probRunList = 4;
zcdo = zcomp( probSetPrm, algoSetPrm, default_solverPrm, zcompPrm );
%
zrun__finish; % Does nothing?
%
msg( __FILE__, __LINE__, "Goodbye." );
printf("\n\n");
