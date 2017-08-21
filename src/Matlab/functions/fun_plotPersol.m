function [  ] = fun_plotPersol( hopfbranch, ind_pext )
%FUN_PLOTPERSOL Summary of this function goes here
%   Detailed explanation goes here

stablty=0;  %unstable 
ind_stabchange=[];
    
for i=1:length(hopfbranch.point)
    if stablty==0 &&(max(hopfbranch.point(i).stability.mu(:))) < 1.01
        stablty=1;   %we accept it as stable
        ind_stabchange=[ind_stabchange i]
    elseif stablty==1 && ((max(hopfbranch.point(i).stability.mu(:))) > 1.01)
        stablty=0;   % if point was stable, but now exceeds limit, we are unstable againg
        ind_stabchange=[ind_stabchange i]
    end    
end

[xm,ym]=df_measr(1,hopfbranch);
ym.subfield='mu';
br_plot(hopfbranch,[],ym,'b');
xlabel('point')
ylabel('Re \lambda')
hopfbranch_psol.point

extval_Hopf=zeros(2,length(hopfbranch.point));
for i=1:length(hopfbranch.point)
   
    extval_Hopf(1,i)=hopfbranch.point(i).parameter(ind_freePam);
    extval_Hopf(2,i)=max(hopfbranch.point(i).profile(2,:) - hopfbranch.point(i).profile(3,:));
    extval_Hopf(3,i)=min(hopfbranch.point(i).profile(2,:) - hopfbranch.point(i).profile(3,:));
end


end



