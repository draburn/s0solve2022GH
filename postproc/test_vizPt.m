clear;
%
test_studyPt;
%close("all");
commondefs;
thisFile = "test_vizPt";
figIndex = 0;
%
tic();
figIndex++; figure(figIndex);
vizPt_vs(studyPtDat,"deltaNorm","omegalin",prm);
figIndex++; figure(figIndex);
vizPt_vs(studyPtDat,"deltaNorm","omega",prm);
prm.figIndex0 = figIndex;
vizPt_curve(studyPtDat,[4,2],prm);
toc;
