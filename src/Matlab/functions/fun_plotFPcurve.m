function [ fix_pot ] = fun_plotFPcurve( fixpointbranch,ind,ind_pext,fignr )
%FUN_PLOTFPCURVE Summary of this function goes here
%   Detailed explanation goes here

%fixpointbranch
fix_pot = zeros(3,length(fixpointbranch.point));
for i=1:length(fixpointbranch.point)
    if max(real(fixpointbranch.point(i).stability.l0)) > 0
        fix_pot(3,i)=1; %branch unstable
    else
        fix_pot(3,i)=-1; %branch stable
    end
    fix_pot(1,i)=fixpointbranch.point(i).parameter(ind_pext);
    fix_pot(2,i)=fixpointbranch.point(i).x(2) - fixpointbranch.point(i).x(3);
end

crossing=[];
for i=1:length(fixpointbranch.point)-1
    if fix_pot(3,i)~=fix_pot(3,i+1)
        crossing = [crossing i];
    end
end

stab=fix_pot(3,1)   %-1 if first point is stable
figure(fignr)
%clf()
markers = [1 crossing length(fixpointbranch.point)];

fireplot=0

if fireplot==0
    
    for i=1:length(markers)-1
        hold on;
        start=markers(i)
        ende=markers(i+1);

        if stab==1
            plot(fix_pot(1,start:ende),fix_pot(2,start:ende)./560,'--m')
        else
            plot(fix_pot(1,start:ende),fix_pot(2,start:ende)./560,'b')
        end
        stab=stab*-1;


    end

elseif fireplot==1
    sig = @(v) 2*2.5 ./ (1+exp(560*6e-3)*exp(-v));
    
    for i=1:length(markers)-1
        hold on;
        start=markers(i)
        ende=markers(i+1);

        if stab==1
            plot(fix_pot(1,start:ende),sig(fix_pot(2,start:ende)),'--m')
        else
            plot(fix_pot(1,start:ende),sig(fix_pot(2,start:ende)),'b')
        end
        stab=stab*-1;
    end
end

 
end

