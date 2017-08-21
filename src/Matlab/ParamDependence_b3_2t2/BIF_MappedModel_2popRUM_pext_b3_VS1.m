
%% pext as primary bifurcation parameter, b3 as secondary bif parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
% ATTENTION: EQUATIONS WERE MODIFIED, SO THAT B1 tunes EIN and external input
%            Also modified was the mapping according to the mapping theory (Version2)

% we gradually increase Nii and motivate that by more and more considering
% collaterals within the inhibitiory faculty


   
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
pext=0;                % signal arriving from other areas
b=taue/taui;            % ration of tau_e to tau_i
b1=0.0;                   % Switch for 2popmodel: if 1, EI is included
b2=1e10;                % not important in this system
b3=0.2;                   % switch for self connectivity of Nii: if 1, selfconn is off
Nep=135;                % Conn gain from py to ei
%Npe=0.8*Nep;            % Conn gain from ei to py
Npe=108
Npi=0.25*Nep;           % Conn gain from ii to py
Nip=0.25*Nep;           % Conn gain from py to ii
Npp=113.4               % self connectivity ei: 
Nii=100;           % self connectivity ii: 
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
% dy(6) = He*taue*r*       (b1*Nep*      sig(x(2)-x(3)) + kp1*n_ext_ei + b1*pext)                                 - 2*x(6)   - x(1);
% dy(7) = He*taue*r*       (b1*Npe*   sig(x(1)     ) + kp3*n_ext_py +(1-b1)*pext + (1-b1)*Npp*sig(x(2)-x(3)) ) - 2*x(7)   - x(2);
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
%     He*taue*r*       (b1*Nep*      sig(x(2,1)-x(3,1)) + kp1*n_ext_ei + b1*pext)                                 - 2*x(6,1)   - x(1,1);...
%     He*taue*r*       (b1*Npe*   sig(x(1,1)     ) + kp3*n_ext_py +(1-b1)*pext + (1-b1)*Npp*sig(x(2,1)-x(3,1)) ) - 2*x(7,1)   - x(2,1);...
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
    par(1)*par(3)*par(19)*    (par(15)*par(7)*      sig(x(2,1)-x(3,1),par) + par(9)*par(12) + par(15)*par(6))                                                - 2*x(6,1)   - x(1,1);...
    par(1)*par(3)*par(19)*    (par(15)*par(8)* sig(x(1,1),par) + par(11)*par(14) + (1-par(15))*par(6) + (1-par(15))*par(23)*sig(x(2,1)-x(3,1),par) ) - 2*x(7,1)   - x(2,1);...
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
freepammin= -180.0;
freepammax = 250.0;
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
freePamRange(2) = 0.50;
freePamRange(3) = 200;
indexHopf       = 188;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_1 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pext,freePamRange,100,0 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Continuation of second Hopf bifurcation in a periodic solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind

freePamRange(1) = -80;
freePamRange(2) = 1.0;
freePamRange(3) = 250;
indexHopf       = 254;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_2 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pext,freePamRange,120,0 );


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
title(['Npp:' num2str(Npp) ' | Nii:' num2str(Nii*(1-b3)) ' | b_1: ' num2str(b1) ' | b_3: ' num2str(b3)])
%axis([-100 250 -0.005 0.03])

print -dpdf 'BIF_map_b3is0k2' -r1000 
bif_map_b3is0k2={fixpointbranch}
save('bif_map_b3is0k2.mat','bif_map_b3is0k2')

% bif_RUM_default={fixpointbranch,hopfbranch_psol_1}
% save('bif_RUM_default.mat','bif_RUM_default')



% or just:
load('bif_RUM_default.mat')
fixpointbranch=bif_RUM_default{1}
hopfbranch_psol_1=bif_RUM_default{2}












%%
% print firing rate fixed point curves for various Nii values
figure(12);clf();
hold on;

load('bif_map_b3is1k0.mat')
fixpointbranch_100=bif_map_b3is1k0{1}
ind=temp_calcind(fixpointbranch_100,ind_pext)
fun_plotFPcurve(fixpointbranch_100,ind,ind_pext,12)
load('bif_map_b3is0k95.mat')
fixpointbranch_95=bif_map_b3is0k95{1}
ind=temp_calcind(fixpointbranch_95,ind_pext)
fun_plotFPcurve(fixpointbranch_95,ind,ind_pext,12)
load('bif_map_b3is0k85.mat')
fixpointbranch_85=bif_map_b3is0k85{1}
ind=temp_calcind(fixpointbranch_85,ind_pext)
fun_plotFPcurve(fixpointbranch_85,ind,ind_pext,12)
load('bif_map_b3is0k8.mat')
fixpointbranch_8=bif_map_b3is0k8{1}
ind=temp_calcind(fixpointbranch_8,ind_pext)
fun_plotFPcurve(fixpointbranch_8,ind,ind_pext,12)
load('bif_map_b3is0k7.mat')
fixpointbranch_7=bif_map_b3is0k7{1}
ind=temp_calcind(fixpointbranch_7,ind_pext)
fun_plotFPcurve(fixpointbranch_7,ind,ind_pext,12)
load('bif_map_b3is0k6.mat')
fixpointbranch_6=bif_map_b3is0k6{1}
ind=temp_calcind(fixpointbranch_6,ind_pext)
fun_plotFPcurve(fixpointbranch_6,ind,ind_pext,12)
load('bif_map_b3is0k5.mat')
fixpointbranch_5=bif_map_b3is0k5{1}
ind=temp_calcind(fixpointbranch_5,ind_pext)
fun_plotFPcurve(fixpointbranch_5,ind,ind_pext,12)
load('bif_map_b3is0k4.mat')
fixpointbranch_4=bif_map_b3is0k4{1}
ind=temp_calcind(fixpointbranch_4,ind_pext)
fun_plotFPcurve(fixpointbranch_4,ind,ind_pext,12)
load('bif_map_b3is0k3.mat')
fixpointbranch_3=bif_map_b3is0k3{1}
ind=temp_calcind(fixpointbranch_3,ind_pext)
fun_plotFPcurve(fixpointbranch_3,ind,ind_pext,12)
load('bif_map_b3is0k2.mat')
fixpointbranch_2=bif_map_b3is0k2{1}
ind=temp_calcind(fixpointbranch_2,ind_pext)
fun_plotFPcurve(fixpointbranch_2,ind,ind_pext,12)


legend('b_3 = 1.0','b_3 = 0.95','b_3 = 0.85','b_3 = 0.8','b_3 = 0.7','b_3 = 0.6','b_3 = 0.5','b_3 = 0.4','b_3 = 0.3','b_3 = 0.2')
% print -dpdf 'Nii_saturation_firerate_He325_Hi22' -r1000 

xlabel('input firing rate')
ylabel('output firing rate')
%ylabel('pyramidal cell membran potential')
title(['Npp:' num2str(Npp) ' | H_e: ' num2str(He*1000) 'mV | H_i: ' num2str(Hi*1000) 'mV'])


%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       b3=0.85
%
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - lower fold - in b3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

% load('bif_map_b3is0k85.mat')
% fixpointbranch=bif_map_b3is0k85{1}
% 

ind_pext
ind_b3=17;

lowerfold_branch_b3=df_brnch(funcs,[ind_pext,ind_b3],'fold'); 

pext_min= -80.0;
pext_max = 400;
pext_step = 1.0;


b3_min = 0.0;
b3_max = 1.0;
b3_step = 0.05;

ind
indexLowerFold       = 323;
lowerfoldpoint=p_tofold(funcs,fixpointbranch.point(indexLowerFold));

lowerfold_branch_b3.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
lowerfold_branch_b3.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
lowerfold_branch_b3.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
lowerfold_branch_b3.point=lowerfoldpoint;
lowerfold_branch_b3.method.continuation.plot=1;

lowerfoldpoint.parameter(ind_b3)=lowerfoldpoint.parameter(ind_b3)-0.001; % perturb hopf point
[lowerfoldpoint,success]=p_correc(funcs,lowerfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
lowerfold_branch_b3.point(2)=lowerfoldpoint;                                 % use as second point of hopf branch:


figure(5); clf;
[lowerfold_branch_b3,s,f,r]=br_contn(funcs,lowerfold_branch_b3,60)           % continue with plotting hopf branch:
lowerfold_branch_b3=br_rvers(lowerfold_branch_b3)                             % reverse Hopf branch
[lowerfold_branch_b3,s,f,r]=br_contn(funcs,lowerfold_branch_b3,50)           % continue in other direction
xlabel('pext');ylabel('b3');

%%
bif_lfb_b3 = zeros(2,length(lowerfold_branch_b3.point));
for i=1:length(lowerfold_branch_b3.point)
    bif_lfb_b3(1,i)=lowerfold_branch_b3.point(i).parameter(ind_pext);
    bif_lfb_b3(2,i)=lowerfold_branch_b3.point(i).parameter(ind_b3);
end

figure(1); clf();
plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:))


figure(16); clf();
plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:),'b')
hold on;
xlabel('pext');ylabel('b3');
title('Variation in b3 - RUM')
legend('lower fold bifurcation branch')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - upper fold - in b3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

% load('bif_map_b3is0k85.mat')
% fixpointbranch=bif_map_b3is0k85{1}
% 


ind_pext
ind_b3=17;

upperfold_branch=df_brnch(funcs,[ind_pext,ind_b3],'fold'); 

pext_min= -250.0;
pext_max = 400;
pext_step = 1.0;


b3_min = 0.0;
b3_max = 1.0;
b3_step = 0.05;

ind
indexUpperFold       = 496;
upperfoldpoint=p_tofold(funcs,fixpointbranch.point(indexUpperFold));

upperfold_branch.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
upperfold_branch.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
upperfold_branch.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
upperfold_branch.point=upperfoldpoint;
upperfold_branch.method.continuation.plot=1;

upperfoldpoint.parameter(ind_b3)=upperfoldpoint.parameter(ind_b3)-0.0001; % perturb hopf point
[upperfoldpoint,success]=p_correc(funcs,upperfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
upperfold_branch.point(2)=upperfoldpoint;                                 % use as second point of hopf branch:


figure(5); clf;
[upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,320)           % continue with plotting hopf branch:
upperfold_branch=br_rvers(upperfold_branch)                             % reverse Hopf branch
[upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,150)           % continue in other direction
xlabel('pext');ylabel('b3');

%%
bif_ufb_b3 = zeros(2,length(upperfold_branch.point));
for i=1:length(upperfold_branch.point)
    bif_ufb_b3(1,i)=upperfold_branch.point(i).parameter(ind_pext);
    bif_ufb_b3(2,i)=upperfold_branch.point(i).parameter(ind_b3);
end

figure(1); clf();
plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:))


figure(16); clf();
plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:),'-.b')
hold on;
plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:),'b')
xlabel('pext');ylabel('b3');
title('Variation in b3 - RUM')
legend('upper fold bifurcation branch','lower fold bifurcation branch')



% print -dpdf 'RUM_2param_2popNii_b3pext_1' -r1000 
% bif_RUM_pext_b3_1={bif_lfb_b3,bif_ufb_b3}
% save('bif_RUM_pext_b3_1.mat','bif_RUM_pext_b3_1')




























%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   b3=1.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - Hopf branch - in b3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams


% load('bif_map_b3is1k0.mat')
% fixpointbranch=bif_map_b3is1k0{1}
% hopfbranch_psol_1=bif_map_b3is1k0{2}




ind_pext
ind_b3=17;

hopfbif_branch_1=df_brnch(funcs,[ind_pext,ind_b3],'hopf'); % use hopf point as first point of hopf branch:

pext_min= -80.0;
pext_max = 250.0;
pext_step = 1.0;


b3_min = 0;
b3_max = 1;
b3_step = 0.1;

ind
indexHopf       = 210;
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_1.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
hopfbif_branch_1.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
hopfbif_branch_1.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
hopfbif_branch_1.point=hopfpoint;
hopfbif_branch_1.method.continuation.plot=1;


hopfpoint.parameter(ind_b3)=hopfpoint.parameter(ind_b3)-0.0001; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_1.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,200)           % continue with plotting hopf branch:
hopfbif_branch_1=br_rvers(hopfbif_branch_1)                             % reverse Hopf branch
[hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,140)           % continue in other direction
xlabel('pext');ylabel('b3');

%%
bif_hopf_b3 = zeros(2,length(hopfbif_branch_1.point));
for i=1:length(hopfbif_branch_1.point)
    bif_hopf_b3(1,i)=hopfbif_branch_1.point(i).parameter(ind_pext);
    bif_hopf_b3(2,i)=hopfbif_branch_1.point(i).parameter(ind_b3);
end

figure(1); clf();
plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:))


n=1
for i=5:7:length(hopfbif_branch_1.point)
    i
    floq(n,1)=hopfbif_branch_1.point(i).parameter(ind_pext); % bif-parameter #1
    floq(n,2)=hopfbif_branch_1.point(i).parameter(ind_b3);   % bif-parameter #2
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
plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:))
xlabel('pext');ylabel('b3');
title('Variation in b3 - RUM')
legend('hopf bifurcation branch 1')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - Hopf branch  2- in b3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams


% load('bif_map_b3is1k0.mat')
% fixpointbranch=bif_map_b3is1k0{1}
% hopfbranch_psol_1=bif_map_b3is1k0{2}




ind_pext
ind_b3=17;

hopfbif_branch_2=df_brnch(funcs,[ind_pext,ind_b3],'hopf'); % use hopf point as first point of hopf branch:

pext_min= -80.0;
pext_max = 850.0;
pext_step = 5.0;


b3_min = 0;
b3_max = 1;
b3_step = 0.1;

ind
indexHopf       = 736;
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_2.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
hopfbif_branch_2.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
hopfbif_branch_2.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
hopfbif_branch_2.point=hopfpoint;
hopfbif_branch_2.method.continuation.plot=1;


hopfpoint.parameter(ind_b3)=hopfpoint.parameter(ind_b3)-0.0001; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_2.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_2,s,f,r]=br_contn(funcs,hopfbif_branch_2,200)           % continue with plotting hopf branch:
hopfbif_branch_2=br_rvers(hopfbif_branch_2)                             % reverse Hopf branch
[hopfbif_branch_2,s,f,r]=br_contn(funcs,hopfbif_branch_2,140)           % continue in other direction
xlabel('pext');ylabel('b3');

%%
bif_hopf_b3_2 = zeros(2,length(hopfbif_branch_2.point));
for i=1:length(hopfbif_branch_2.point)
    bif_hopf_b3_2(1,i)=hopfbif_branch_2.point(i).parameter(ind_pext);
    bif_hopf_b3_2(2,i)=hopfbif_branch_2.point(i).parameter(ind_b3);
end

figure(1); clf();
plot(bif_hopf_b3_2(1,:), bif_hopf_b3_2(2,:))



n=1
for i=65:4:length(hopfbif_branch_2.point)
    i
    floq(n,1)=hopfbif_branch_2.point(i).parameter(ind_pext); % bif-parameter #1
    floq(n,2)=hopfbif_branch_2.point(i).parameter(ind_b3);   % bif-parameter #2
    hopfpoint=hopfbif_branch_2.point(i);
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
plot(bif_hopf_b3_2(1,:), bif_hopf_b3_2(2,:))
hold on;
plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:))
plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:),'-.b')
plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:),'b')
xlabel('pext');ylabel('b3');
title('Variation in b3 - RUM')
legend('hopf bifurcation branch 1','hopf bifurcation branch 2','upper fold bifurcation branch','lower fold bifurcation branch')



% print -dpdf 'RUM_2param_2popNii_b3pext_allStabilities' -r1000 
% bif_RUM_pext_b3_1={bif_hopf_b3_2,bif_hopf_b3,bif_lfb_b3,bif_ufb_b3}
% save('bif_RUM_pext_b3_1.mat','bif_RUM_pext_b3_1')



load('bif_RUM_pext_b3_1.mat')
bif_hopf_b3_2=bif_RUM_pext_b3_1{1}
bif_hopf_b3=bif_RUM_pext_b3_1{2}
bif_lfb_b3=bif_RUM_pext_b3_1{3}
bif_ufb_b3=bif_RUM_pext_b3_1{4}









% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   2 parameter bifurcation - lower fold - in b3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% ind_pext
% ind_b3=17;
% 
% lowerfold_branch_b3=df_brnch(funcs,[ind_pext,ind_b3],'fold'); 
% 
% pext_min= -80.0;
% pext_max = 200;
% pext_step = 1.0;
% 
% 
% b3_min = 0.0;
% b3_max = 1.0;
% b3_step = 0.05;
% 
% ind
% indexLowerFold       = 272;
% lowerfoldpoint=p_tofold(funcs,fixpointbranch.point(indexLowerFold));
% 
% lowerfold_branch_b3.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
% lowerfold_branch_b3.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
% lowerfold_branch_b3.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
% lowerfold_branch_b3.point=lowerfoldpoint;
% lowerfold_branch_b3.method.continuation.plot=1;
% 
% lowerfoldpoint.parameter(ind_b3)=lowerfoldpoint.parameter(ind_b3)-0.001; % perturb hopf point
% [lowerfoldpoint,success]=p_correc(funcs,lowerfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
% lowerfold_branch_b3.point(2)=lowerfoldpoint;                                 % use as second point of hopf branch:
% 
% 
% figure(5); clf;
% [lowerfold_branch_b3,s,f,r]=br_contn(funcs,lowerfold_branch_b3,100)           % continue with plotting hopf branch:
% lowerfold_branch_b3=br_rvers(lowerfold_branch_b3)                             % reverse Hopf branch
% [lowerfold_branch_b3,s,f,r]=br_contn(funcs,lowerfold_branch_b3,100)           % continue in other direction
% xlabel('pext');ylabel('b3');
% 
% 
% bif_lfb_b3 = zeros(2,length(lowerfold_branch_b3.point));
% for i=1:length(lowerfold_branch_b3.point)
%     bif_lfb_b3(1,i)=lowerfold_branch_b3.point(i).parameter(ind_pext);
%     bif_lfb_b3(2,i)=lowerfold_branch_b3.point(i).parameter(ind_b3);
% end
% 
% figure(1); clf();
% plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:))
% 
% 
% figure(16); clf();
% plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:))
% hold on;
% plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:))
% hold on;
% axis([-100 180 -0.1 1.1])
% %plot([-100 20],[1.0 1.0],'k.-')
% xlabel('pext');ylabel('b3');
% title('Variation in b3 - RUM')
% legend('lower fold bifurcation branch','hopf bifurcation branch')
% 
% 
% % print -dpdf 'RUM_2PamBIF_pext_b3_NppNii0k3' -r1000 
% % bif_RUM_pext_b3_NppNii0k3={bif_lfb_b3,bif_hopf_b3}
% % save('bif_RUM_pext_b3_NppNii0k3.mat','bif_RUM_pext_b3_NppNii0k3')
% % 
% 
% lowerfold_branch_b3=br_stabl(funcs,lowerfold_branch_b3,0,0);
% figure(); clf;
% [xm,ym]=df_measr(1,lowerfold_branch_b3); % plot stability versus point number:
% ym.subfield='l0';
% br_plot(lowerfold_branch_b3,[],ym,'c');
% ym.subfield='l1';
% br_plot(lowerfold_branch_b3,[],ym,'b');
% xlabel('point number along branch');ylabel('\Re(\lambda)');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   2 parameter bifurcation - upper fold - in b3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% ind_pext
% ind_b3=17;
% 
% upperfold_branch=df_brnch(funcs,[ind_pext,ind_b3],'fold'); 
% 
% pext_min= -80.0;
% pext_max = 200;
% pext_step = 1.0;
% 
% 
% b3_min = 0.0;
% b3_max = 1.0;
% b3_step = 0.05;
% 
% ind
% indexUpperFold       = 393;
% upperfoldpoint=p_tofold(funcs,fixpointbranch.point(indexUpperFold));
% 
% upperfold_branch.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
% upperfold_branch.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
% upperfold_branch.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
% upperfold_branch.point=upperfoldpoint;
% upperfold_branch.method.continuation.plot=1;
% 
% upperfoldpoint.parameter(ind_b3)=upperfoldpoint.parameter(ind_b3)-0.0001; % perturb hopf point
% [upperfoldpoint,success]=p_correc(funcs,upperfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
% upperfold_branch.point(2)=upperfoldpoint;                                 % use as second point of hopf branch:
% 
% 
% figure(5); clf;
% [upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,220)           % continue with plotting hopf branch:
% upperfold_branch=br_rvers(upperfold_branch)                             % reverse Hopf branch
% [upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,450)           % continue in other direction
% xlabel('pext');ylabel('b3');
% 
% %%
% bif_ufb_b3 = zeros(2,length(upperfold_branch.point));
% for i=1:length(upperfold_branch.point)
%     bif_ufb_b3(1,i)=upperfold_branch.point(i).parameter(ind_pext);
%     bif_ufb_b3(2,i)=upperfold_branch.point(i).parameter(ind_b3);
% end
% 
% figure(1); clf();
% plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:))
% 
% 
% figure(16); clf();
% plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:),'-.b')
% hold on;
% plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:))
% plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:),'g')
% hold on;
% axis([-100 180 -0.1 1.3])
% %plot([-100 20],[1.0 1.0],'k.-')
% xlabel('pext');ylabel('b3');
% title('Variation in b3 - RUM')
% legend('upper fold bifurcation branch','lower fold bifurcation branch','hopf bifurcation branch')
% 
% 
% % print -dpdf 'RUM_2PamBIF_pext_VS2_Npp125' -r1000 
% % bif_RUM_pext_b3_VS2_Npp125={bif_ufb_b3,bif_lfb_b3,bif_hopf_b3}
% % save('bif_RUM_pext_b3_VS2_Npp125.mat','bif_RUM_pext_b3_VS2_Npp125')
% % 
% 
% load('bif_RUM_pext_b3_VS2_Npp125.mat','bif_RUM_pext_b3_VS2_Npp125')
% bif_ufb_b3=bif_RUM_pext_b3_VS2_Npp125{1}
% bif_lfb_b3=bif_RUM_pext_b3_VS2_Npp125{2}
% bif_hopf_b3=bif_RUM_pext_b3_VS2_Npp125{3}
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %       b=0
% %
% %%   2 parameter bifurcation - Hopf branch - in b3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% 
% 
% load('bif_RUM_b30_Npp125.mat')
% fixpointbranch=bif_RUM_b30_Npp125{1}
% hopfbranch_psol_1=bif_RUM_b30_Npp125{2}
% 
% 
% 
% 
% %%
% 
% ind_pext
% ind_b3=17;
% 
% hopfbif_branch_1=df_brnch(funcs,[ind_pext,ind_b3],'hopf'); % use hopf point as first point of hopf branch:
% 
% pext_min= -80.0;
% pext_max = 500.0;
% pext_step = 0.5;
% 
% 
% b3_min = 0;
% b3_max = 1;
% b3_step = 0.1;
% 
% ind
% indexHopf       = 299;
% hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));
% 
% hopfbif_branch_1.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
% hopfbif_branch_1.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
% hopfbif_branch_1.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
% hopfbif_branch_1.point=hopfpoint;
% hopfbif_branch_1.method.continuation.plot=1;
% 
% 
% hopfpoint.parameter(ind_b3)=hopfpoint.parameter(ind_b3)+0.0001; % perturb hopf point
% [hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
% hopfbif_branch_1.point(2)=hopfpoint;                                 % use as second point of hopf branch:
% 
% 
% figure(6); clf;
% [hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,510)           % continue with plotting hopf branch:
% hopfbif_branch_1=br_rvers(hopfbif_branch_1)                             % reverse Hopf branch
% [hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,140)           % continue in other direction
% xlabel('pext');ylabel('b3');
% 
% %%
% bif_hopf_b3_0 = zeros(2,length(hopfbif_branch_1.point));
% for i=1:length(hopfbif_branch_1.point)
%     bif_hopf_b3_0(1,i)=hopfbif_branch_1.point(i).parameter(ind_pext);
%     bif_hopf_b3_0(2,i)=hopfbif_branch_1.point(i).parameter(ind_b3);
% end
% 
% figure(1); clf();
% plot(bif_hopf_b3_0(1,:), bif_hopf_b3_0(2,:))
% 
% 
% figure(16); clf();
% 
% plot(bif_hopf_b3_0(1,:), bif_hopf_b3_0(2,:),'m')
% hold on;
% plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:),'-.b')
% plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:))
% plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:),'g')
% 
% axis([-100 220 -0.1 1.1])
% xlabel('pext');ylabel('b3');
% title('Variation in b3 - RUM')
% legend('upper fold bifurcation branch','lower fold bifurcation branch','hopf bifurcation branch_b3_1','hopf bifurcation branch_b3_0_1')
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %       b=0
% %
% %%   2 parameter bifurcation - 2nd Hopf branch - in b3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% 
% 
% load('bif_RUM_b30_Npp125.mat')
% fixpointbranch=bif_RUM_b30_Npp125{1}
% hopfbranch_psol_1=bif_RUM_b30_Npp125{2}
% %%
% 
% 
% 
% 
% 
% ind_pext
% ind_b3=17;
% 
% hopfbif_branch_b0_2=df_brnch(funcs,[ind_pext,ind_b3],'hopf'); % use hopf point as first point of hopf branch:
% 
% pext_min= -80.0;
% pext_max = 12000.0;
% pext_step = 50.0;
% 
% 
% b3_min = 0;
% b3_max = 1;
% b3_step = 0.1;
% 
% ind
% indexHopf       = 790;
% hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));
% 
% hopfbif_branch_b0_2.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
% hopfbif_branch_b0_2.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
% hopfbif_branch_b0_2.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
% hopfbif_branch_b0_2.point=hopfpoint;
% hopfbif_branch_b0_2.method.continuation.plot=1;
% 
% 
% hopfpoint.parameter(ind_b3)=hopfpoint.parameter(ind_b3)+0.0001; % perturb hopf point
% [hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
% hopfbif_branch_b0_2.point(2)=hopfpoint;                                 % use as second point of hopf branch:
% 
% 
% figure(6); clf;
% [hopfbif_branch_b0_2,s,f,r]=br_contn(funcs,hopfbif_branch_b0_2,510)           % continue with plotting hopf branch:
% hopfbif_branch_b0_2=br_rvers(hopfbif_branch_b0_2)                             % reverse Hopf branch
% [hopfbif_branch_b0_2,s,f,r]=br_contn(funcs,hopfbif_branch_b0_2,140)           % continue in other direction
% xlabel('pext');ylabel('b3');
% %%
% 
% bif_hopf_b3_0_2 = zeros(2,length(hopfbif_branch_b0_2.point));
% for i=1:length(hopfbif_branch_b0_2.point)
%     bif_hopf_b3_0_2(1,i)=hopfbif_branch_b0_2.point(i).parameter(ind_pext);
%     bif_hopf_b3_0_2(2,i)=hopfbif_branch_b0_2.point(i).parameter(ind_b3);
% end
% 
% figure(1); clf();
% plot(bif_hopf_b3_0_2(1,:), bif_hopf_b3_0_2(2,:))
% 
% 
% figure(16); clf();
% plot(bif_hopf_b3_0_2(1,:), bif_hopf_b3_0_2(2,:),'k')
% hold on;
% plot(bif_hopf_b3_0(1,:), bif_hopf_b3_0(2,:),'m')
% plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:),'-.b')
% plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:))
% plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:),'g')
% 
% axis([-100 12000 -0.1 1.1])
% %plot([-100 20],[1.0 1.0],'k.-')
% xlabel('pext');ylabel('b3');
% title('Variation in b3 - RUM')
% legend('upper fold bifurcation branch','lower fold bifurcation branch','hopf bifurcation branch_b3_1','hopf bifurcation branch_b3_0_1','hopf bifurcation branch_b3_0_2')
% 
% 
% % print -dpdf 'RUM_2PamBIF_pext_VS2_Npp125_longhopf' -r1000 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   2 parameter bifurcation - hopf-fold bif - in b3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% load('bif_RUM_b30_Npp125.mat')
% fixpointbranch=bif_RUM_b30_Npp125{1}
% hopfbranch_psol_1=bif_RUM_b30_Npp125{2}
% 
% 
% ind_pext
% ind_b3=17;
% 
% 
% 
% indexhopffoldbif       = 36; % get parameter as minimal val (fold) from hopfbranch_psol_1 (generate hopf fold for He=3.4208mV)
% %hopfbranch_psol_1.point(36).parameter(6)
% %initialize branch
% 
% pext_max = 150;
% pext_step = 0.1;
% 
% 
% hopfbranch_psol_1.parameter.max_step=[ind_pext,pext_max]; % remove step size restriction
% [pfuncs,pbranch,suc]=SetupPOfold(funcs,hopfbranch_psol_1,indexhopffoldbif,...
%     'contpar',[ind_pext,ind_b3],'dir',ind_pext,'step',pext_step,'print_residual_info',1);
% if suc
%     disp('POFold initialization finished');
% else
%     warning('POFold initialization failed');
% end
% 
% figure(2);clf()
% pbranch=br_contn(pfuncs,pbranch,100);
%    
% 
% bif_hfb_b3 = zeros(2,length(pbranch.point));
% for i=1:length(pbranch.point)
%     bif_hfb_b3(1,i)=pbranch.point(i).parameter(ind_pext);
%     bif_hfb_b3(2,i)=pbranch.point(i).parameter(ind_b3);
% end
% 
% figure(3);clf();
% plot(bif_hfb_b3(1,:),bif_hfb_b3(2,:),'b')
% 
% 
% 
% figure(16); clf();
% plot(bif_hfb_b3(1,:),bif_hfb_b3(2,:),'k')
% hold on;
% plot(bif_hopf_b3_0(1,:), bif_hopf_b3_0(2,:),'m')
% plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:),'-.b')
% plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:))
% plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:),'g')
% 
% axis([-100 220 -0.1 1.1])
% %plot([-100 20],[1.0 1.0],'k.-')
% xlabel('pext');ylabel('b3');
% title('Variation in b3 - RUM')
% legend('hopf fold','hopf bifurcation branch b3=0 1','upper fold bifurcation branch','lower fold bifurcation branch','hopf bifurcation branch b3=1')
% 
% 
% % print -dpdf 'RUM_2PamBIF_pext_VS2_Npp125' -r1000 
% % bif_RUM_pext_b3_VS2_Npp125={bif_ufb_b3,bif_lfb_b3,bif_hopf_b3,bif_hopf_b3_0,bif_hopf_b3_0_2,bif_hfb_b3}
% % save('bif_RUM_pext_b3_VS2_Npp125.mat','bif_RUM_pext_b3_VS2_Npp125')
% % 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% %
% %% b3=0.6
% %
% %
% %
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   2 parameter bifurcation - lower fold - in b3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% load('bif_RUM_b30k6_Npp125.mat')
% fixpointbranch=bif_RUM_b30k6_Npp125{1}
% 
% 
% 
% 
% 
% 
% ind_pext
% ind_b3=17;
% 
% lowerfold_branch_b3=df_brnch(funcs,[ind_pext,ind_b3],'fold'); 
% 
% pext_min= -80.0;
% pext_max = 200;
% pext_step = 1.0;
% 
% 
% b3_min = 0.0;
% b3_max = 1.0;
% b3_step = 0.05;
% 
% ind
% indexLowerFold       = 340;
% lowerfoldpoint=p_tofold(funcs,fixpointbranch.point(indexLowerFold));
% 
% lowerfold_branch_b3.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
% lowerfold_branch_b3.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
% lowerfold_branch_b3.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
% lowerfold_branch_b3.point=lowerfoldpoint;
% lowerfold_branch_b3.method.continuation.plot=1;
% 
% lowerfoldpoint.parameter(ind_b3)=lowerfoldpoint.parameter(ind_b3)+0.001; % perturb hopf point
% [lowerfoldpoint,success]=p_correc(funcs,lowerfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
% lowerfold_branch_b3.point(2)=lowerfoldpoint;                                 % use as second point of hopf branch:
% 
% 
% figure(5); clf;
% [lowerfold_branch_b3,s,f,r]=br_contn(funcs,lowerfold_branch_b3,100)           % continue with plotting hopf branch:
% lowerfold_branch_b3=br_rvers(lowerfold_branch_b3)                             % reverse Hopf branch
% [lowerfold_branch_b3,s,f,r]=br_contn(funcs,lowerfold_branch_b3,50)           % continue in other direction
% xlabel('pext');ylabel('b3');
% 
% %%
% bif_lfb_b3_0k6 = zeros(2,length(lowerfold_branch_b3.point));
% for i=1:length(lowerfold_branch_b3.point)
%     bif_lfb_b3_0k6(1,i)=lowerfold_branch_b3.point(i).parameter(ind_pext);
%     bif_lfb_b3_0k6(2,i)=lowerfold_branch_b3.point(i).parameter(ind_b3);
% end
% 
% figure(1); clf();
% plot(bif_lfb_b3_0k6(1,:), bif_lfb_b3_0k6(2,:),'r')
% 
% 
% 
% figure(16); clf();
% plot(bif_lfb_b3_0k6(1,:), bif_lfb_b3_0k6(2,:),'r')
% %plot(bif_hfb_b3(1,:),bif_hfb_b3(2,:),'k')
% hold on;
% plot(bif_hopf_b3_0(1,:), bif_hopf_b3_0(2,:),'m')
% %plot(bif_hopf_b3_0_2(1,:), bif_hopf_b3_0_2(2,:),'m')
% plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:),'-.b')
% plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:))
% plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:),'g')
% 
% axis([-100 220 -0.1 1.1])
% %plot([-100 20],[1.0 1.0],'k.-')
% xlabel('pext');ylabel('b3');
% title('Variation in b3 - RUM')
% legend('lower fold bifurcation branch b3=0.6','hopf bifurcation branch b3=0 #1','upper fold bifurcation branch','lower fold bifurcation branch','hopf bifurcation branch b3=1')
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   2 parameter bifurcation - upper fold - in b3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% load('bif_RUM_b30k6_Npp125.mat')
% fixpointbranch=bif_RUM_b30k6_Npp125{1}
% 
% 
% 
% ind_pext
% ind_b3=17;
% 
% upperfold_branch=df_brnch(funcs,[ind_pext,ind_b3],'fold'); 
% 
% pext_min= -80.0;
% pext_max = 200;
% pext_step = 1.0;
% 
% 
% b3_min = 0.0;
% b3_max = 1.0;
% b3_step = 0.05;
% 
% ind
% indexUpperFold       = 361;
% upperfoldpoint=p_tofold(funcs,fixpointbranch.point(indexUpperFold));
% 
% upperfold_branch.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_b3 b3_min]']';
% upperfold_branch.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_b3 b3_max]']';
% upperfold_branch.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_b3 b3_step]']';
% upperfold_branch.point=upperfoldpoint;
% upperfold_branch.method.continuation.plot=1;
% 
% upperfoldpoint.parameter(ind_b3)=upperfoldpoint.parameter(ind_b3)-0.0001; % perturb hopf point
% [upperfoldpoint,success]=p_correc(funcs,upperfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
% upperfold_branch.point(2)=upperfoldpoint;                                 % use as second point of hopf branch:
% 
% 
% figure(5); clf;
% [upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,40)           % continue with plotting hopf branch:
% upperfold_branch=br_rvers(upperfold_branch)                             % reverse Hopf branch
% [upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,200)           % continue in other direction
% xlabel('pext');ylabel('b3');
% 
% %%
% bif_ufb_b3_0k6 = zeros(2,length(upperfold_branch.point));
% for i=1:length(upperfold_branch.point)
%     bif_ufb_b3_0k6(1,i)=upperfold_branch.point(i).parameter(ind_pext);
%     bif_ufb_b3_0k6(2,i)=upperfold_branch.point(i).parameter(ind_b3);
% end
% 
% figure(1); clf();
% plot(bif_ufb_b3_0k6(1,:), bif_ufb_b3_0k6(2,:))
% 
% figure(16); clf();
% plot(bif_ufb_b3_0k6(1,:), bif_ufb_b3_0k6(2,:),'-.r')
% hold on;
% plot(bif_lfb_b3_0k6(1,:), bif_lfb_b3_0k6(2,:),'r')
% %plot(bif_hfb_b3(1,:),bif_hfb_b3(2,:),'k')
% plot(bif_hopf_b3_0(1,:), bif_hopf_b3_0(2,:),'m')
% %plot(bif_hopf_b3_0_2(1,:), bif_hopf_b3_0_2(2,:),'m')
% plot(bif_ufb_b3(1,:), bif_ufb_b3(2,:),'-.b')
% plot(bif_lfb_b3(1,:), bif_lfb_b3(2,:))
% plot(bif_hopf_b3(1,:), bif_hopf_b3(2,:),'g')
% 
% axis([-100 220 -0.1 1.1])
% %plot([-100 20],[1.0 1.0],'k.-')
% xlabel('pext');ylabel('b3');
% title('Variation in b3 - RUM')
% legend('upper fold bifurcation branch b3=0.6','lower fold bifurcation branch b3=0.6','hopf bifurcation branch b3=0 #1','upper fold bifurcation branch','lower fold bifurcation branch','hopf bifurcation branch b3=1')
% 
% % bif_RUM_pext_b3_VS2_Npp125={bif_ufb_b3_0k6,bif_lfb_b3_0k6,bif_ufb_b3,bif_lfb_b3,bif_hopf_b3,bif_hopf_b3_0,bif_hopf_b3_0_2}
% % save('bif_RUM_pext_b3_VS2_Npp125.mat','bif_RUM_pext_b3_VS2_Npp125')
% 
% 
% % print -dpdf 'RUM_2PamBIF_pext_VS2_Npp125' -r1000 
% % 
% % load('bif_RUM_pext_b3_VS2_Npp125.mat','bif_RUM_pext_b3_VS2_Npp125')
% % bif_ufb_b3_0k6=bif_RUM_pext_b3_VS2_Npp125{1}
% % bif_lfb_b3_0k6=bif_RUM_pext_b3_VS2_Npp125{2}
% % bif_ufb_b3=bif_RUM_pext_b3_VS2_Npp125{3}
% % bif_lfb_b3=bif_RUM_pext_b3_VS2_Npp125{4}
% % bif_hopf_b3=bif_RUM_pext_b3_VS2_Npp125{5}
% % bif_hopf_b3_0=bif_RUM_pext_b3_VS2_Npp125{6}
% % bif_hopf_b3_0_2=bif_RUM_pext_b3_VS2_Npp125{7}
