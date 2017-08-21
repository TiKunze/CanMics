function [ floquet_mp ] = fun_determ_Floquet_HpfBifBranch( funcs,hopfpoint,fixpointbranch,indexHopf,index_freePam,freePamRange,NumPoints)

%   This file sumes up the necessary steps to continue a peridioc solution
%   for conitnuation of a Hopf cycle and returning the Floquet multiplier
%   to determine the type of hopf bifurcation

flag_newhheur =1;

%hopf=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopf=hopfpoint;%p_tohopf(funcs,fixpointbranch.point(indexHopf));

% correct the point to find the actual hopf bifurcation
method=df_mthod(funcs,'hopf',flag_newhheur); % get hopf calculation method parameters:
method.stability.minimal_real_part=-100;
[hopf,success]=p_correc(funcs,hopf,index_freePam,[],method.point) % correct hopf

if ~success
    [hopf,success]=p_correc(funcs,hopf,index_freePam,[],method.point)
end


if ~success
    error('Hopf correction not successfull')
end


% convert hopf point to the nearby limit cycle psol
amp=1e-2; degree=3; intervals=18;
[psol_point,stepcond]=p_topsol(funcs,hopf,amp,degree,intervals);  %convert hopf point to the nearby limit cycle psol

% correct the psol
method=df_mthod(funcs,'psol');
[psol_point,success]=p_correc(funcs,psol_point,index_freePam,stepcond,method.point); %freePam (here pext)
if ~success,
    [psol_point,success]=p_correc(funcs,psol_point,index_freePam,stepcond,method.point); %freePam (here pext)
end

if ~success,
    [psol_point,success]=p_correc(funcs,psol_point,index_freePam,stepcond,method.point); %freePam (here pext)
end

if ~success,
    fprint(1,'correction failed\n')
end

% define a psol branch
hopfbranch_psol=df_brnch(funcs,index_freePam,'psol');
freepammin= freePamRange(1)
freepammax = freePamRange(3)
freepamstep = freePamRange(2)

hopfbranch_psol.parameter.min_bound(1,:)=[index_freePam freepammin];
hopfbranch_psol.parameter.max_bound(1,:)=[index_freePam freepammax];
hopfbranch_psol.parameter.max_step(1,:)=[index_freePam freepamstep];

hopfbranch_psol.method.continuation.plot=0;
% first point is degenerate 0 amplitude solution
[deg_psol_point,stepcond]=p_topsol(funcs,hopf,0,degree,intervals);
hopfbranch_psol.point=deg_psol_point;
hopfbranch_psol.point(2)=psol_point;



% continue in one direction:
[hopfbranch_psol,s,f,r]=br_contn(funcs,hopfbranch_psol,NumPoints)



% calculate stability
hopfbranch_psol=br_stabl(funcs,hopfbranch_psol,0,1);


floquet_mp=hopfbranch_psol.point(end).stability.mu;


end

