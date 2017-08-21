function [ index, points,dirs ] = fun_detectzerocross( stream, branchpoints, indices )
%DETECTZEROCROSS Summary of this function goes here
%   this script comes in handy when you want to detect the zero crossings
%   of eigenvalues. The crossings are detected at the parameter stream and
%   as ponit indices in the branch

index=[];
points=[];
dirs=[];

for i=1:length(stream)-1
   if (stream(i)<0 && stream(i+1)>0) 
       dirs=[dirs +1];
       index=[index indices(i)];
       points=[points branchpoints(i)];
   elseif (stream(i)>0 && stream(i+1)<0)
       dirs=[dirs -1];
       index=[index indices(i)];
       points=[points branchpoints(i)];    
   end
end


end

