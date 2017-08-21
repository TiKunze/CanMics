function [ cutted_branch ] = fun_cutbranch( branch,begin,ende )
%FUN_CUTPATH Summary of this function goes here
%   Detailed explanation goes here

cutted_branch.method=branch.method;
cutted_branch.parameter=branch.parameter;
cutted_branch.point=branch.point(begin:ende);

end

