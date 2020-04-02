clear;
%
test_studyPt;
%close("all");
commondefs;
thisFile = "test_vizPt";
%
tic();
prm.figIndex0 = 0;
figure;
vizPt_vs(studyPtDat,"deltaNorm","omega",prm);
figure;
vizPt_curve1d(studyPtDat,1,prm);
figure;
vizPt_curve1d(studyPtDat,8,prm);
figure;
vizPt_curve1d(studyPtDat,[1,8],prm);
prm.figIndex0=4;
vizPt_curve(studyPtDat,[1,8],prm);
toc;
