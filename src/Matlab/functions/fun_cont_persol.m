function [ hopfbranch_psol ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,index_freePam,freePamRange,NumPoints,plotyn )
%CONT_HOPFBR Summary of this function goes here
%   This file sumes up the necessary steps to continue a peridioc solution
%   for conitnuation of a HOpf cycle

% it is not generic, but needs to be adapted if used

flag_newhheur =1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Continuation of first Hopf bifurcation in a periodic solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hopf=p_tohopf(funcs,fixpointbranch.point(indexHopf));

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

% plot the eigenvalues to verify
hopf.stability=p_stabil(funcs,ststpoint,method.stability); % compute stability of hopf point
%figure(); clf;
%p_splot(hopf);                     % plot stability of hopf point


% convert hopf point to the nearby limit cycle psol

ststpoint=fixpointbranch.point(500)
amp=1e-2; degree=3; intervals=18;
[psol_point,stepcond]=p_topsol(funcs,ststpoint,amp,degree,intervals);  %convert hopf point to the nearby limit cycle psol

% correct the psol
method=df_mthod(funcs,'psol');
[psol_point,success]=p_correc(funcs,psol_point,index_freePam,stepcond,method.point) %freePam (here pext)
if ~success,
    fprint(1,'correction failed\n')
end

% define a psol branch
hopfbranch_psol=df_brnch(funcs,index_freePam,'psol')
freepammin= freePamRange(1);
freepammax = freePamRange(3);
freepamstep = freePamRange(2);

hopfbranch_psol.parameter.min_bound(1,:)=[index_freePam freepammin];
hopfbranch_psol.parameter.max_bound(1,:)=[index_freePam freepammax];
hopfbranch_psol.parameter.max_step(1,:)=[index_freePam freepamstep];
hopfbranch_psol.method.continuation.plot=0;
if plotyn==1
    hopfbranch_psol.method.continuation.plot=1;
end

% first point is degenerate 0 amplitude solution
[deg_psol_point,stepcond]=p_topsol(funcs,hopf,0,degree,intervals);
hopfbranch_psol.point=deg_psol_point;
hopfbranch_psol.point(2)=psol_point;



% continue in one direction:
if hopfbranch_psol.method.continuation.plot == 1
    figure();clf();
end
[hopfbranch_psol,s,f,r]=br_contn(funcs,hopfbranch_psol,NumPoints)



% calculate stability
hopfbranch_psol=br_stabl(funcs,hopfbranch_psol,0,1);


% plot magnitude of floquet exponents along the branch
figure(100); clf();
subplot(221)

[xm,ym]=df_measr(1,hopfbranch_psol);
ym.subfield='mu';
br_plot(hopfbranch_psol,[],ym,'b');
xlabel('point')
ylabel('Re \lambda')


% plot frequency along the branch
subplot(222)
freqplot_br(hopfbranch_psol) 


% plot point to parameter mapping
subplot(223)
for i=1:length(hopfbranch_psol.point)
    pam(i)=hopfbranch_psol.point(i).parameter(index_freePam);
end
plot(1:length(hopfbranch_psol.point),pam)
xlabel('point number')
ylabel('parameter')
title('point - parameter mapping')
end

