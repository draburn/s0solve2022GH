clear;
%
test_studyPt;
%close("all");
commondefs;
thisFile = "test_vizPt";
%
tic();
figure(1);
vizPt_vs(studyPtDat,"deltaNorm","omega",prm);
figure(2);
vizPt_curve1d(studyPtDat,1,prm);
figure(3);
vizPt_curve1d(studyPtDat,8,prm);
figure(4);
vizPt_curve1d(studyPtDat,[1,8],prm);
prm.figIndex0 = 4;
vizPt_curve(studyPtDat,[1,8],prm);
toc;
