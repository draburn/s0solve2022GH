thisFile = "groot1d__finish";
%
%
msg_main( verbLev, thisFile, __LINE__, sprintf( ...
  "Final:  %4d  %6.2f  %8.2e  %s", ...
   fevalCount, ...
   time()-startTime, ...
   fNormMin, ...
   retcode2str(retCode) )  );
%
%
datOut.funchF = funchF;
datOut.x1 = x1;
datOut.x2 = x2;
datOut.prm = prm;
datOut.startTime = startTime;
% Do NOT (blindly) copy datIn! Required memory could get too large!
%
datOut.verbLev = verbLev;
datOut.reportInterval = reportInterval;
%
datOut.fNormTol = fNormTol;
datOut.exeTimeLimit = exeTimeLimit;
datOut.fevalCountLimit = fevalCountLimit;
%
datOut.fevalCount = fevalCount;
datOut.fNormMin = fNormMin;
%
datOut.xVals = xVals;
datOut.fVals = fVals;
datOut.evalIndex_sorted = evalIndex_sorted;
datOut.xVals_sorted = xVals_sorted;
datOut.fVals_sorted = fVals_sorted;
%
%
return;
