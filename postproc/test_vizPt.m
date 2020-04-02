clear;
%
test_studyPt;
close("all");
commondefs;
thisFile = "test_vizPt";
%
tic();
prm.figIndex0 = 0;
figure(1);
vizPt_vs(studyPtDat,"deltaNorm","omega",prm);
figure(2);
vizPt_curve(studyPtDat,1,prm);
prm.figIndex0=2;
vizPt_curve(studyPtDat,[1,2],prm);
toc;
