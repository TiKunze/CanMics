
%% pext as primary bifurcation parameter, b1 as secondary bif parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Current version: in version 3 we double scaled Npe,Nep,Npp with b1: once
% in the preparation and then in the equations again. in this version I
% directly put the definitions for Nep,Npe,and Npp into the right side of
% the equation. Note, however, that our mapping theory still assumes an
% alpha=(1-b1)
%
%
% ATTENTION: EQUATIONS WERE MODIFIED, SO THAT B1 tunes EIN and external input
%            Also modified was the mapping according to the text (Version2)

%Np'e=Npe*(1-alpha)/(1+alpha*(Me/Mp))
%Np'p'=Npe*alpha/(1+alpha*(Me/Mp)) + Nep*alpha/(alpha+Mp/Me)
%Nep'=Nep

%Np'e=Npe*(b1)/(1+(1-b1)*(Me/Mp))
%Np'p'=Npe*(1-b1)/(1+(1-b1)*(Me/Mp)) + Nep*(1-b1)/((1-b1)+Mp/Me)
%Nep'=Nep



   
%            b1 stays b1   
%            b2 gets  b1
%            Nep is multiplied by b1
%            b3 at Npp gets b1


run init_ddebiftool;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of user functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
He=3.25e-3;             % V
Hi=22e-3;               % V
taue=0.01;              % s
taui=0.02;              % s
kp1=0;                  % connectivity of external noise to excitatory cells
kp2=0;                  % connectivity of external noise to inhibitory cells
kp3=0;                  % connectivity of external noise to pyramidal cells
n_ext_ii=0;             % external noise to inhibitory cells
n_ext_ei=0;             % external noise to excitatory cells
n_ext_py=0;             % external noise to pyramidal cells
pext=50;                % signal arriving from other areas
b=taue/taui;            % ration of tau_e to tau_i
b1=0.10;                % 3pop switch: if 1, EI is included and receives input
b2=1e10;                % not important in this system
b3=1;                   % switch for self connectivity of Nii: if 1, selfconn is off
Nep=135;                % Conn gain from py to ei

Npe=0.8*Nep;           % Conn gain from ei to py 
Npi=0.25*Nep;           % Conn gain from ii to py
Nip=0.25*Nep;           % Conn gain from py to ii
Npp=0*1e5;                % DUMMY: HERE:directly defined in equations   self connectivity ei: 
Nii=0*100000.95*Npi;           % self connectivity ii: assumption Nii~Npi
e0=2.5;                 % 1/s
r=560;                  % 1/V
v0=6e-3;                % V

par = zeros(24,1);

par(1) = He;
par(2) = Hi;
par(3) = taue;
par(4) = taui;
par(5) = b;
par(6) = pext;
par(7) = Nep;
par(8) = Npe;
par(9) = kp1;            % ext to EI
par(10) = kp2;            % ext to II
par(11) = kp3;            % ext to Py
par(12) = n_ext_ei;       % ext to EI
par(13) = n_ext_ii;       % ext to II
par(14) = n_ext_py;       % ext to Py
par(15) = b1;           % Switch for EI: if 1, EI is included
par(16) = b2;           % switch for input: if 0, input is fed to Py
par(17) = b3;           % switch for self connectivity: if 1, selfconn is off
par(18) = e0;
par(19) = r;
par(20) = v0;
par(21) = Npi;
par(22) = Nip;
par(23) = Npp;
par(24) = Nii;
par(25) = 0.0001; 



% %generic mapped equation system:
% dy(1) = x(6);
% dy(2) = x(7);
% dy(3) = x(8);
% dy(4) = x(9);1
% dy(5) = x(10);
% dy(6) = He*taue*r*       (Nep'*      sig(x(2)-x(3)) + kp1*n_ext_ei + b1*pext)                                 - 2*x(6)   - x(1);
% dy(7) = He*taue*r*       (Np'e*   sig(x(1)     ) + kp3*n_ext_py +(1-b1)*pext + Np'p'*sig(x(2)-x(3)) ) - 2*x(7)   - x(2);
% dy(8) = Hi*taui*b^2*r*Npi*            sig(x(4)-x(5))                                                           - 2*b*x(8) - b^2*x(3);
% dy(9) = He*taue*r*(Nip*sig(x(2)-x(3)) + kp2*n_ext_ii )                                           - 2*x(9)   - x(4);
% dy(10)= Hi*taui*b^2*r*Nii*(1-b3)*sig(x(4)-x(5))                                              - 2*b*x(10) - b^2*x(5);
% 
% MappedModelGeneric_sys_rhs=@(x,par)[...
%     x(6,1);...
%     x(7,1);...
%     x(8,1);...
%     x(9,1);...
%     x(10,1);...
%     He*taue*r*       (Nep'*      sig(x(2,1)-x(3,1)) + kp1*n_ext_ei + b1*pext)                                 - 2*x(6,1)   - x(1,1);...
%     He*taue*r*       (Np'e*   sig(x(1,1)     ) + kp3*n_ext_py +(1-b1)*pext + Np'p'*sig(x(2,1)-x(3,1)) ) - 2*x(7,1)   - x(2,1);...
%     Hi*taui*b^2*r*Npi*            sig(x(4,1)-x(5,1))                                                           - 2*b*x(8,1) - b^2*x(3,1);...
%     He*taue*r*(Nip*sig(x(2,1)-x(3,1)) + kp2*n_ext_ii )                                           - 2*x(9,1)   - x(4,1);...
%     Hi*taui*b^2*r*Nii*(1-b3)*sig(x(4,1)-x(5,1))                                              - 2*b*x(10,1) - b^2*x(5,1)]

sig = @(v,par) 2*par(18) ./ (1+exp(par(19)*par(20))*exp(-v));

MappedModelGeneric_sys_rhs=@(x,par)[...
    x(6,1);...
    x(7,1);...
    x(8,1);...
    x(9,1);...
    x(10,1);...
    par(1)*par(3)*par(19)*    (par(7)*      sig(x(2,1)-x(3,1),par) + par(9)*par(12) + par(15)*par(6))                                                - 2*x(6,1)   - x(1,1);...
    par(1)*par(3)*par(19)*    (par(8)*(par(15))/(1+(1-par(15))*(0.25)) * sig(x(1,1),par) + par(11)*par(14) + (1-par(15))*par(6) + (par(8)*(1-par(15))/(1+(1-par(15))*(0.25)) + par(7)*(1-par(15))/((1-par(15))+4)) *sig(x(2,1)-x(3,1),par) ) - 2*x(7,1)   - x(2,1);...
    par(2)*par(4)*par(5)^2*   par(19)*par(21)*            sig(x(4,1)-x(5,1),par)                                                                     - 2*par(5)*x(8,1) - par(5)^2*x(3,1);...
    par(1)*par(3)*par(19)*    (par(22)*sig(x(2,1)-x(3,1),par) + par(10)*par(13))                                                                     - 2*x(9,1)   - x(4,1);...
    par(2)*par(4)*par(5)^2*par(19)*par(24)*(1-par(17))    *sig(x(4,1)-x(5,1),par)                                                                    - 2*par(5)*x(10,1) - par(5)^2*x(5,1)];




% Delays
neuron_tau=@()[25];

% Bifurcation parameter
ind_pext=6;

funcs=set_funcs(...
    'sys_rhs',MappedModelGeneric_sys_rhs,...
    'sys_tau',neuron_tau)

%% Initialise the first guessed fixed point
stst.kind='stst';
stst.parameter=par';
stst.x=zeros(10,1);

% Calculate Fixpoint curve 

flag_newhheur=1; % use the new steplength heuristic (Verheyden_2007)
method=df_mthod(funcs,'stst',flag_newhheur);
method.stability.minimal_real_part=-40
[stst,success] = p_correc(funcs,stst,[],[],method.point)
stst.x
if ~success,
    fprintf(1,'correction failed\n')
end

% compute its stability:
stst.stability = p_stabil(funcs,stst,method.stability)
method.stability.minimal_real_part=-100; 
stst.stability=p_stabil(funcs,stst,method.stability); % recompute stability:
figure(1); clf;
p_splot(stst); 




% Initialize branch of trivial equilibria
% get an empty branch with ind_pext as a free parameter:
fixpointbranch=df_brnch(funcs,ind_pext,'stst')

% set bounds for continuation parameter
freepammin= -50.0;
freepammax = 850.0;
freepamstep = 0.50;

fixpointbranch.parameter.min_bound(1,:)=[ind_pext freepammin];
fixpointbranch.parameter.max_bound(1,:)=[ind_pext freepammax];
fixpointbranch.parameter.max_step(1,:)=[ind_pext freepamstep];
% use stst as a first branch point:
fixpointbranch.point=stst;


%  Extend and continue branch of trivial equilibria
stst.parameter(ind_pext)=stst.parameter(ind_pext)+0.0150;
method=df_mthod(funcs,'stst')
[stst,success]=p_correc(funcs,stst,[],[],method.point)

% use as a second branch point:
fixpointbranch.point(2)=stst;
fixpointbranch.method.continuation.plot=0;     %switch off plotting
fixpointbranch.method.stability.minimal_real_part=-100;


figure(1);clf();
% continue in one direction:
[fixpointbranch,s,f,r]=br_contn(funcs,fixpointbranch,3300)
% turn the branch around:
fixpointbranch=br_rvers(fixpointbranch);
% continue in the other direction:
[fixpointbranch,s,f,r]=br_contn(funcs,fixpointbranch,3000)
fixpointbranch=br_rvers(fixpointbranch);


% compute stability of the branch
% the first 0 is how many points to skip between stab calculations
% the second 0 is to not recalculate any stability already present
fixpointbranch=br_stabl(funcs,fixpointbranch,0,0);


% obtain suitable scalar measures to plot stability along branch:
[xm,ym]=df_measr(1,fixpointbranch)
figure(1); clf;
br_plot(fixpointbranch,xm,ym,'b');     % plot stability along branch:
ym.subfield='l0';
br_plot(fixpointbranch,xm,ym,'c');     % l0: rightmost roots
plot([freepammin freepammax],[0 0],'-.');
xlabel('free parameter');ylabel('\Re\lambda');

%
%Plot Eigenvalues for specific point
ind_point=31;
figure(1);clf();
p_splot(fixpointbranch.point(ind_point));
title(['point: ' num2str(ind_point) ' | parameter:' num2str(fixpointbranch.point(ind_point).parameter(ind_pext))])
%

% plot stability versus point number:
figure(1); clf;
br_plot(fixpointbranch,[],ym,'b');
br_plot(fixpointbranch,[],ym,'b.');
plot([0 30],[0 0],'-.');
xlabel('point number along branch');ylabel('\Re(\lambda)');

% Plot all Eigenvalues along branch and along free parameter

branch_summar=zeros(22,length(fixpointbranch.point));
for i=1:length(fixpointbranch.point)
    branch_summar(1,i)=fixpointbranch.point(i).parameter(ind_pext);   % parameter value
    branch_summar(2,i)=i;                                         % point number
    branch_summar(3:12,i)=fixpointbranch.point(i).x;                     % fix point vals
    branch_summar(13:22,i)=fixpointbranch.point(i).stability.l0;         % stability of eigen values
end

%Detect zero crossings
ind=[[],[],[]];

for i=1:10
    [ind_temp,points_temp, dirs ]=fun_detectzerocross(real(branch_summar(12+i,:)),branch_summar(2,:),branch_summar(1,:));
    ind=[ind [dirs;ind_temp;points_temp]];
    %points=[points points_temp];
end

figure(1);clf();
for i=1:10
    subplot(2,10,i)
    plot(branch_summar(1,:),real(branch_summar(12+i,:))); 
    hold on;
    plot([branch_summar(1,1) branch_summar(1,end)], [0 0],'k.-');
    subplot(2,10,10+i)
    xlabel('input firing rate')
    ylabel('\Re(\lambda)')
    plot(branch_summar(2,:),real(branch_summar(12+i,:)));
    hold on;
    plot([branch_summar(1,1) branch_summar(1,end)], [0 0],'k.-');
    xlabel('point number')
    ylabel('\Re(\lambda)')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Continuation of first Hopf bifurcation in a periodic solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind

freePamRange(1) = -50;
freePamRange(2) = 7.5;
freePamRange(3) = 750;
indexHopf       = 399;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_1 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pext,freePamRange,300,1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Continuation of second Hopf bifurcation in a periodic solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind

freePamRange(1) = 10;
freePamRange(2) = 5.0;
freePamRange(3) = 850;
indexHopf       = 1812;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_2 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pext,freePamRange,250,1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Put Branches together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first hopf branch

[ extval_Hopf_1 ] = fun_prepareBranchPlot( hopfbranch_psol_1 , ind_pext);
[ extval_Hopf_2 ] = fun_prepareBranchPlot( hopfbranch_psol_2 , ind_pext);

figure(11);clf();
%plot(fix_pot(1,:),fix_pot(2,:)./560,'b')
fun_plotFPcurve(fixpointbranch,ind,ind_pext,11)

hold on;
plot(extval_Hopf_1(1,:),extval_Hopf_1(2,:)./560,'r');
plot(extval_Hopf_1(1,:),extval_Hopf_1(3,:)./560,'r');

plot(extval_Hopf_2(1,:),extval_Hopf_2(2,:)./560,'g');
plot(extval_Hopf_2(1,:),extval_Hopf_2(3,:)./560,'g');

xlabel('input firing rate')
ylabel('pyrapot')
Npe=par(8)*(par(15))/(1+(1-par(15))*(0.25))
Npp=(par(8)*(1-par(15))/(1+(1-par(15))*(0.25)) + par(7)*(1-par(15))/((1-par(15))+4))
title(['Npp:' num2str(Npp) ' | Npe:' num2str(Npe) ' | Nep:' num2str(Nep) ' | b_1: ' num2str(b1)])

print -dpdf 'BIF_map_b1is0k1' -r1000 
bif_map_b0k1={fixpointbranch,hopfbranch_psol_1}
save('bif_map_b0k1.mat','bif_map_b0k1')

% bif_RUM_default={fixpointbranch,hopfbranch_psol_1}
% save('bif_RUM_default.mat','bif_RUM_default')



% or just:
load('bif_RUM_default.mat')
fixpointbranch=bif_RUM_default{1}
hopfbranch_psol_1=bif_RUM_default{2}
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - Hopf branch - in b1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

% load('bif_map_b1k0.mat')
% fixpointbranch=bif_map_b1k0{1}
% hopfbranch_psol_1=bif_map_b1k0{2}

ind_pext
ind_b1=15;

hopfbif_branch_1=df_brnch(funcs,[ind_pext,ind_b1],'hopf'); % use hopf point as first point of hopf branch:

pext_min= -80.0;
pext_max = 800.0;
pext_step = 3.0;


b1_min = 0;
b1_max = 1;
b1_step = 0.1;

ind
indexHopf       = 326;
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_1.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b1 b1_min]']';
hopfbif_branch_1.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b1 b1_max]']';
hopfbif_branch_1.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b1 b1_step]']';
hopfbif_branch_1.point=hopfpoint;
hopfbif_branch_1.method.continuation.plot=1;


hopfpoint.parameter(ind_b1)=hopfpoint.parameter(ind_b1)-0.0001; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_1.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,500)           % continue with plotting hopf branch:
hopfbif_branch_1=br_rvers(hopfbif_branch_1)                             % reverse Hopf branch
[hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,140)           % continue in other direction
xlabel('pext');ylabel('b1');

%%
bif_hopf_b1 = zeros(2,length(hopfbif_branch_1.point));
for i=1:length(hopfbif_branch_1.point)
    bif_hopf_b1(1,i)=hopfbif_branch_1.point(i).parameter(ind_pext);
    bif_hopf_b1(2,i)=hopfbif_branch_1.point(i).parameter(ind_b1);
end

figure(1); clf();
plot(bif_hopf_b1(1,:), bif_hopf_b1(2,:))




n=1
for i=10:10:length(hopfbif_branch_1.point)
    i
    floq(n,1)=hopfbif_branch_1.point(i).parameter(ind_pext);
    floq(n,2)=hopfbif_branch_1.point(i).parameter(ind_b1);
    hopfpoint=hopfbif_branch_1.point(i);
    PamRange(1)=hopfpoint.parameter(ind_pext) - 0.2*abs(hopfpoint.parameter(ind_pext));
    PamRange(3)=hopfpoint.parameter(ind_pext) + 0.2*abs(hopfpoint.parameter(ind_pext));
    PamRange(2)=(abs(PamRange(3))-abs(PamRange(2)))/30;
    temp=fun_determ_Floquet_HpfBifBranch( funcs,hopfpoint,fixpointbranch,indexHopf,ind_pext,PamRange,20)
    if length(temp)==2
        floq(n,3:4)=temp;
    elseif length(temp)==4
        floq(n,3:6)=temp;
        display "complex"
    else
        display "AAAAAAAAAAAAAAAAAAAAAAA"
    end
    n=n+1
    i
end

figure(16)
clf()
hold on;
for i=1:length(floq(:,1))
    if max(floq(i,3:4))>1.005
        plot(floq(i,1),floq(i,2),'.r','MarkerSize',20)
    elseif max(floq(i,3:4))<1.001    
        plot(floq(i,1),floq(i,2),'.g','MarkerSize',20)
    else
        plot(floq(i,1),floq(i,2),'.m','MarkerSize',20)
    end
end








figure(16); clf();
plot(bif_hopf_b1(1,:), bif_hopf_b1(2,:))
hold on;
axis([-100 150 -0.1 1.1])
%plot([-100 20],[1.0 1.0],'k.-')
xlabel('pext');ylabel('b1');
title('Variation in b1 - RUM')
legend('hopf bifurcation branch')













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - lower fold - in b1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams


% load('bif_map_b1k0.mat')
% fixpointbranch=bif_map_b1k0{1}
% hopfbranch_psol_1=bif_map_b1k0{2}

ind_pext
ind_b1=15;

lowerfold_branch_b1=df_brnch(funcs,[ind_pext,ind_b1],'fold'); 

pext_min= -80.0;
pext_max = 800;
pext_step = 1.0;


b1_min = 0.0;
b1_max = 1.0;
b1_step = 0.5;

ind
indexLowerFold       = 172;
lowerfoldpoint=p_tofold(funcs,fixpointbranch.point(indexLowerFold));

lowerfold_branch_b1.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b1 b1_min]']';
lowerfold_branch_b1.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b1 b1_max]']';
lowerfold_branch_b1.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b1 b1_step]']';
lowerfold_branch_b1.point=lowerfoldpoint;
lowerfold_branch_b1.method.continuation.plot=1;

lowerfoldpoint.parameter(ind_b1)=lowerfoldpoint.parameter(ind_b1)-0.001; % perturb hopf point
[lowerfoldpoint,success]=p_correc(funcs,lowerfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
lowerfold_branch_b1.point(2)=lowerfoldpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[lowerfold_branch_b1,s,f,r]=br_contn(funcs,lowerfold_branch_b1,100)           % continue with plotting hopf branch:
lowerfold_branch_b1=br_rvers(lowerfold_branch_b1)                             % reverse Hopf branch
[lowerfold_branch_b1,s,f,r]=br_contn(funcs,lowerfold_branch_b1,100)           % continue in other direction
xlabel('pext');ylabel('b1');

%%
bif_lfb_b1 = zeros(2,length(lowerfold_branch_b1.point));
for i=1:length(lowerfold_branch_b1.point)
    bif_lfb_b1(1,i)=lowerfold_branch_b1.point(i).parameter(ind_pext);
    bif_lfb_b1(2,i)=lowerfold_branch_b1.point(i).parameter(ind_b1);
end

figure(1); clf();
plot(bif_lfb_b1(1,:), bif_lfb_b1(2,:))





figure(16); clf();
plot(bif_lfb_b1(1,:), bif_lfb_b1(2,:),'b')
hold on;
plot(bif_hopf_b1(1,:), bif_hopf_b1(2,:))
hold on;
axis([-100 180 -0.1 1.1])
%plot([-100 20],[1.0 1.0],'k.-')
xlabel('pext');ylabel('b1');
title('Variation in b1 - RUM')
legend('lower fold bifurcation branch','hopf bifurcation branch')


% print -dpdf 'RUM_2PamBIF_pext_b1_NppNii0k3' -r1000 
% bif_RUM_pext_b1_NppNii0k3={bif_lfb_b1,bif_hopf_b1}
% save('bif_RUM_pext_b1_NppNii0k3.mat','bif_RUM_pext_b1_NppNii0k3')
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - upper fold - in b1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams



% load('bif_map_b1k0.mat')
% fixpointbranch=bif_map_b1k0{1}
% hopfbranch_psol_1=bif_map_b1k0{2}


ind_pext
ind_b1=15;

upperfold_branch=df_brnch(funcs,[ind_pext,ind_b1],'fold'); 

pext_min= -80.0;
pext_max = 800;
pext_step = 1.0;


b1_min = 0.0;
b1_max = 1.0;
b1_step = 0.05;

ind
indexUpperFold       = 293;
upperfoldpoint=p_tofold(funcs,fixpointbranch.point(indexUpperFold));

upperfold_branch.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b1 b1_min]']';
upperfold_branch.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b1 b1_max]']';
upperfold_branch.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b1 b1_step]']';
upperfold_branch.point=upperfoldpoint;
upperfold_branch.method.continuation.plot=1;

upperfoldpoint.parameter(ind_b1)=upperfoldpoint.parameter(ind_b1)-0.0001; % perturb hopf point
[upperfoldpoint,success]=p_correc(funcs,upperfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
upperfold_branch.point(2)=upperfoldpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,220)           % continue with plotting hopf branch:
upperfold_branch=br_rvers(upperfold_branch)                             % reverse Hopf branch
[upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,450)           % continue in other direction
xlabel('pext');ylabel('b1');

%%
bif_ufb_b1 = zeros(2,length(upperfold_branch.point));
for i=1:length(upperfold_branch.point)
    bif_ufb_b1(1,i)=upperfold_branch.point(i).parameter(ind_pext);
    bif_ufb_b1(2,i)=upperfold_branch.point(i).parameter(ind_b1);
end

figure(1); clf();
plot(bif_ufb_b1(1,:), bif_ufb_b1(2,:))


figure(16); clf();
plot(bif_ufb_b1(1,:), bif_ufb_b1(2,:),'-.b')
hold on;
plot(bif_lfb_b1(1,:), bif_lfb_b1(2,:),'b')
plot(bif_hopf_b1(1,:), bif_hopf_b1(2,:),'g')
hold on;
axis([-100 180 -0.1 1.3])
%plot([-100 20],[1.0 1.0],'k.-')
xlabel('pext');ylabel('b1');
title('Variation in b1 - RUM')
legend('upper fold bifurcation branch','lower fold bifurcation branch','hopf bifurcation branch')


% print -dpdf 'RUM_2PamBIF_pext_VS4_Npp113k4' -r1000 
% bif_RUM_pext_b1_VS4_Npp113k4={bif_ufb_b1,bif_lfb_b1,bif_hopf_b1,floq}
% save('bif_RUM_pext_b1_VS4_Npp113k4.mat','bif_RUM_pext_b1_VS4_Npp113k4')
% % 
% 
% load('bif_RUM_pext_b1_VS2_Npp125.mat','bif_RUM_pext_b1_VS2_Npp125')
% bif_ufb_b1=bif_RUM_pext_b1_VS2_Npp125{1}
% bif_lfb_b1=bif_RUM_pext_b1_VS2_Npp125{2}
% bif_hopf_b1=bif_RUM_pext_b1_VS2_Npp125{3}


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       b=0
%
%%   2 parameter bifurcation - Hopf branch - in b1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams


% 
% load('bif_map_b0k0.mat')
% fixpointbranch=bif_map_b0k0{1}
% hopfbranch_psol_1=bif_map_b0k0{2}



%%

ind_pext
ind_b1=15;

hopfbif_branch_1=df_brnch(funcs,[ind_pext,ind_b1],'hopf'); % use hopf point as first point of hopf branch:

pext_min= -80.0;
pext_max = 800.0;
pext_step = 0.5;


b1_min = 0;
b1_max = 1;
b1_step = 0.1;

ind
indexHopf       = 210;
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_1.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b1 b1_min]']';
hopfbif_branch_1.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b1 b1_max]']';
hopfbif_branch_1.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b1 b1_step]']';
hopfbif_branch_1.point=hopfpoint;
hopfbif_branch_1.method.continuation.plot=1;


hopfpoint.parameter(ind_b1)=hopfpoint.parameter(ind_b1)+0.0001; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_1.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,510)           % continue with plotting hopf branch:
hopfbif_branch_1=br_rvers(hopfbif_branch_1)                             % reverse Hopf branch
[hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,140)           % continue in other direction
xlabel('pext');ylabel('b1');

%%
bif_hopf_b1_0 = zeros(2,length(hopfbif_branch_1.point));
for i=1:length(hopfbif_branch_1.point)
    bif_hopf_b1_0(1,i)=hopfbif_branch_1.point(i).parameter(ind_pext);
    bif_hopf_b1_0(2,i)=hopfbif_branch_1.point(i).parameter(ind_b1);
end

figure(1); clf();
plot(bif_hopf_b1_0(1,:), bif_hopf_b1_0(2,:))



n=1
for i=5:10:length(hopfbif_branch_1.point)
    i
    floq(n,1)=hopfbif_branch_1.point(i).parameter(ind_pext);
    floq(n,2)=hopfbif_branch_1.point(i).parameter(ind_b1);
    hopfpoint=hopfbif_branch_1.point(i);
    PamRange(1)=hopfpoint.parameter(ind_pext) - 0.2*abs(hopfpoint.parameter(ind_pext));
    PamRange(3)=hopfpoint.parameter(ind_pext) + 0.2*abs(hopfpoint.parameter(ind_pext));
    PamRange(2)=(abs(PamRange(3))-abs(PamRange(2)))/30;
    temp=fun_determ_Floquet_HpfBifBranch( funcs,hopfpoint,fixpointbranch,indexHopf,ind_pext,PamRange,20)
    if length(temp)==2
        floq(n,3:4)=temp;
    elseif length(temp)==4
        floq(n,3:6)=temp;
        display "complex"
    else
        display "AAAAAAAAAAAAAAAAAAAAAAA"
    end
    n=n+1
    i
end

figure(16)
clf()
hold on;
for i=1:length(floq(:,1))
    if max(floq(i,3:4))>1.005
        plot(floq(i,1),floq(i,2),'.r','MarkerSize',20)
    elseif max(floq(i,3:4))<1.001    
        plot(floq(i,1),floq(i,2),'.g','MarkerSize',20)
    else
        plot(floq(i,1),floq(i,2),'.m','MarkerSize',20)
    end
end




figure(16); clf();

plot(bif_hopf_b1_0(1,:), bif_hopf_b1_0(2,:),'m')
hold on;
plot(bif_ufb_b1(1,:), bif_ufb_b1(2,:),'-.b')
plot(bif_lfb_b1(1,:), bif_lfb_b1(2,:))
plot(bif_hopf_b1(1,:), bif_hopf_b1(2,:),'g')

axis([-100 220 -0.1 1.1])
xlabel('pext');ylabel('b1');
title('Variation in b1 - RUM')
legend('upper fold bifurcation branch','lower fold bifurcation branch','hopf bifurcation branch_b1_1','hopf bifurcation branch_b1_0_1')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       b=0
%
%%   2 parameter bifurcation - 2nd Hopf branch - in b1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams



% 
% load('bif_map_b0k0.mat')
% fixpointbranch=bif_map_b0k0{1}
% hopfbranch_psol_1=bif_map_b0k0{2}
%%





ind_pext
ind_b1=15;

hopfbif_branch_b0_2=df_brnch(funcs,[ind_pext,ind_b1],'hopf'); % use hopf point as first point of hopf branch:

pext_min= -80.0;
pext_max = 750.0;
pext_step = 50.0;


b1_min = 0;
b1_max = 1;
b1_step = 0.1;

ind
indexHopf       = 736;
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_b0_2.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b1 b1_min]']';
hopfbif_branch_b0_2.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b1 b1_max]']';
hopfbif_branch_b0_2.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b1 b1_step]']';
hopfbif_branch_b0_2.point=hopfpoint;
hopfbif_branch_b0_2.method.continuation.plot=1;


hopfpoint.parameter(ind_b1)=hopfpoint.parameter(ind_b1)+0.0001; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_b0_2.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_b0_2,s,f,r]=br_contn(funcs,hopfbif_branch_b0_2,510)           % continue with plotting hopf branch:
hopfbif_branch_b0_2=br_rvers(hopfbif_branch_b0_2)                             % reverse Hopf branch
[hopfbif_branch_b0_2,s,f,r]=br_contn(funcs,hopfbif_branch_b0_2,140)           % continue in other direction
xlabel('pext');ylabel('b1');
%%

bif_hopf_b1_0_2 = zeros(2,length(hopfbif_branch_b0_2.point));
for i=1:length(hopfbif_branch_b0_2.point)
    bif_hopf_b1_0_2(1,i)=hopfbif_branch_b0_2.point(i).parameter(ind_pext);
    bif_hopf_b1_0_2(2,i)=hopfbif_branch_b0_2.point(i).parameter(ind_b1);
end

figure(1); clf();
plot(bif_hopf_b1_0_2(1,:), bif_hopf_b1_0_2(2,:))


n=1
for i=5:4:length(hopfbif_branch_b0_2.point)
    i
    floq(n,1)=hopfbif_branch_b0_2.point(i).parameter(ind_pext);
    floq(n,2)=hopfbif_branch_b0_2.point(i).parameter(ind_b1);
    hopfpoint=hopfbif_branch_b0_2.point(i);
    PamRange(1)=hopfpoint.parameter(ind_pext) - 0.2*abs(hopfpoint.parameter(ind_pext));
    PamRange(3)=hopfpoint.parameter(ind_pext) + 0.2*abs(hopfpoint.parameter(ind_pext));
    PamRange(2)=(abs(PamRange(3))-abs(PamRange(2)))/30;
    temp=fun_determ_Floquet_HpfBifBranch( funcs,hopfpoint,fixpointbranch,indexHopf,ind_pext,PamRange,20)
    if length(temp)==2
        floq(n,3:4)=temp;
    elseif length(temp)==4
        floq(n,3:6)=temp;
        display "complex"
    else
        display "AAAAAAAAAAAAAAAAAAAAAAA"
    end
    n=n+1
    i
end

figure(16)
clf()
hold on;
for i=1:length(floq(:,1))
    if max(floq(i,3:4))>1.005
        plot(floq(i,1),floq(i,2),'.r','MarkerSize',20)
    elseif max(floq(i,3:4))<1.001    
        plot(floq(i,1),floq(i,2),'.g','MarkerSize',20)
    else
        plot(floq(i,1),floq(i,2),'.m','MarkerSize',20)
    end
end




figure(16); clf();
plot(bif_hopf_b1_0_2(1,:), bif_hopf_b1_0_2(2,:),'k')
hold on;
plot(bif_hopf_b1_0(1,:), bif_hopf_b1_0(2,:),'m')
plot(bif_ufb_b1(1,:), bif_ufb_b1(2,:),'-.b')
plot(bif_lfb_b1(1,:), bif_lfb_b1(2,:),'b')
plot(bif_hopf_b1(1,:), bif_hopf_b1(2,:),'g')

axis([-100 750 -0.1 1.1])
%plot([-100 20],[1.0 1.0],'k.-')
xlabel('pext');ylabel('b1');
title('Variation in b1 - RUM')
legend('hopf bifurcation branch b1 0_2','hopf bifurcation branch b1 0_1','upper fold bifurcation branch','lower fold bifurcation branch','hopf bifurcation branch b1 1')

print -dpdf 'RUM_2PamBIF_pext_VS4_Npp113k4_allstabilities' -r1000 
bif_RUM_pext_b1_VS4_Npp113k4={bif_hopf_b1_0_2,bif_hopf_b1_0,bif_ufb_b1,bif_lfb_b1,bif_hopf_b1,floq}
save('bif_RUM_pext_b1_VS4_Npp113k4.mat','bif_RUM_pext_b1_VS4_Npp113k4')
 









