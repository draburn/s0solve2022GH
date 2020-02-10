thisFile = "linsolf__finish";
%
%
datOut.funchMatAProd = funchMatAProd;
datOut.vecB = vecB;
datOut.prm = prm;
datOut.thisFile = thisFile;
datOut.startTime = startTime;
%
datOut.verbLev = verbLev;
datOut.reportInterval = reportInterval;
%
datOut.sizeF = sizeF;
datOut.sizeX = sizeX;
datOut.normB = normB;
%
datOut.fracResTol = fracResTol;
datOut.exeTimeLimit = exeTimeLimit;
datOut.numIterLimit = numIterLimit;
%
datOut.gsThresh0 = gsThresh0;
datOut.gsThresh1 = gsThresh1;
%
datOut.numIter = numIter;
datOut.matU = matU;
datOut.matV = matV;
datOut.matW = matW;
datOut.vecX = vecX;
datOut.fracRes = fracRes;
%
%
msg_main( verbLev, thisFile, __LINE__, sprintf( ...
  "Final:  %4d  %6.2f  %8.2e  %s", ...
   numIter, ...
   time()-startTime, ...
   fracRes, ...
   retcode2str(retCode) )  );
%
return;
