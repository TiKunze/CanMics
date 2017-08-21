function [ extval_Hopf ] = fun_prepareBranchPlot( hopfbranch,ind_freePam )
%FUN_PREPAREBRANCHPLOT Summary of this function goes here
%   Detailed explanation goes here

extval_Hopf=zeros(2,length(hopfbranch.point));
for i=1:length(hopfbranch.point)
   
    extval_Hopf(1,i)=hopfbranch.point(i).parameter(ind_freePam);
    extval_Hopf(2,i)=max(hopfbranch.point(i).profile(2,:) - hopfbranch.point(i).profile(3,:));
    extval_Hopf(3,i)=min(hopfbranch.point(i).profile(2,:) - hopfbranch.point(i).profile(3,:));
end


end

