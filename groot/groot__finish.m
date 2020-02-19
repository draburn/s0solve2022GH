thisFile = "groot__finish";
%
%
msg_main( verbLev, thisFile, __LINE__, sprintf( ...
  "Final:  %4d  %6.2f  %8.2e  %s", ...
   numIter, ...
   time()-startTime, ...
   omega, ...
   retcode2str(retCode) )  );
%
%
datOut.vecX0 = vecX0;
datOut.funchF = funchF;
datOut.prm = prm;
datOut.startTime = startTime;
% Do NOT (blindly) copy datIn! Required memory could get too large!
%
datOut.verbLev = verbLev;
datOut.reportInterval = reportInterval;
%
datOut.sizeX = sizeX;
datOut.sizeF = sizeF;
%
datOut.omegaTol = omegaTol;
datOut.exeTimeLimit = exeTimeLimit;
datOut.numIterLimit = numIterLimit;
%
datOut.numIter = numIter;
datOut.omega = omega;
%
datOut.vecX = vecX;
%
%
return;
