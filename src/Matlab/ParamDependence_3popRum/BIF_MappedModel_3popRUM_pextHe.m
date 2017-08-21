
%% pext as primary bifurcation parameter, He as secondary bif parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Prepare work space: First Part: pext as primary bifurcation parameter

run init_ddebiftool;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of user functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
He=3.25e-3;             % V
Hi=22e-3;               % V11.27 - 11.3
taue=0.01;              % s
taui=0.02;              % s
kp1=0;                  % connectivity of external noise to excitatory cells
kp2=0;                  % connectivity of external noise to inhibitory cells
kp3=0;                  % connectivity of external noise to pyramidal cells
n_ext_ii=0;             % external noise to inhibitory cells
n_ext_ei=0;             % external noise to excitatory cells
n_ext_py=0;             % external noise to pyramidal cells
pext=-50;                % signal arriving from other areas
b=taue/taui;            % ration of tau_e to tau_i
b1=1;                   % Switch for EI: if 1, EI is included
b2=1;                   % switch for input: if 0, input is fed to Py
b3=1;                   % switch for self connectivity: if 1, selfconn is off
Nep=135;                % Conn gain from py to ei
Npe=0.8*Nep;            % Conn gain from ei to py
Npi=0.25*Nep;           % Conn gain from ii to py
Nip=0.25*Nep;           % Conn gain from py to ii
Npp=1e6*Nip;           % self connectivity ei: assumption Npp~Nip
Nii=1e6*Npi;           % self connectivity ii: assumption Nii~Npi
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
% dy(2) = x(7);4
% dy(3) = x(8);
% dy(4) = x(9);
% dy(5) = x(10);
% dy(6) = He*taue*r*       (Nep*      sig(x(2)-x(3)) + kp1*n_ext_ei + b2*pext)                                 - 2*x(6)   - x(1);
% dy(7) = He*taue*r*       (b1*Npe*   sig(x(1)     ) + kp3*n_ext_py +(1-b2)*pext + (1-b3)*Npp*sig(x(2)-x(3)) ) - 2*x(7)   - x(2);
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
%     He*taue*r*       (Nep*      sig(x(2,1)-x(3,1)) + kp1*n_ext_ei + b2*pext)                                 - 2*x(6,1)   - x(1,1);...
%     He*taue*r*       (b1*Npe*   sig(x(1,1)     ) + kp3*n_ext_py +(1-b2)*pext + (1-b3)*Npp*sig(x(2,1)-x(3,1)) ) - 2*x(7,1)   - x(2,1);...
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
    par(1)*par(3)*par(19)*    (par(7)*      sig(x(2,1)-x(3,1),par) + par(9)*par(12) + par(16)*par(6))                                                - 2*x(6,1)   - x(1,1);...
    par(1)*par(3)*par(19)*    (par(15)*par(8)* sig(x(1,1),par) + par(11)*par(14) + (1-par(16))*par(6) + (1-par(17))*par(23)*sig(x(2,1)-x(3,1),par) ) - 2*x(7,1)   - x(2,1);...
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

% Initialise the first guessed fixed point
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
freepammin= -250.0;
freepammax = 300.0;
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


%Plot Eigenvalues for specific point
ind_point=222;
figure(1);clf();
p_splot(fixpointbranch.point(ind_point));
title(['point: ' num2str(ind_point) ' | parameter:' num2str(fixpointbranch.point(ind_point).parameter(ind_pext))])


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

freePamRange(1) = -250;
freePamRange(2) = 1.0;
freePamRange(3) = 200;
indexHopf       = 958%1196;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_1 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pext,freePamRange,90,1 );

%hopfbranch_1 = fun_cutbranch( hopfbranch_psol_1,1,84 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Continuation of second Hopf bifurcation in a periodic solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind

freePamRange(1) = -200;
freePamRange(2) = .5;
freePamRange(3) = 250;
indexHopf       = 1614;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_2 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pext,freePamRange,150,1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Put Branches together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first hopf branch
[ extval_Hopf_1 ] = fun_prepareBranchPlot( hopfbranch_psol_1 , ind_pext);
[ extval_Hopf_2 ] = fun_prepareBranchPlot( hopfbranch_psol_2 , ind_pext);


figure(12);clf();
%plot(fix_pot(1,:),fix_pot(2,:)./560,'b')
fun_plotFPcurve(fixpointbranch,ind,ind_pext,12)
hold on;
%fun_plotPersol(hopfbranch_1, ind_pext)

plot(extval_Hopf_1(1,:),extval_Hopf_1(2,:)./560,'r');
plot(extval_Hopf_1(1,:),extval_Hopf_1(3,:)./560,'r');

plot(extval_Hopf_2(1,:),extval_Hopf_2(2,:)./560,'g');
plot(extval_Hopf_2(1,:),extval_Hopf_2(3,:)./560,'g');

xlabel('input firing rate')
ylabel('pyrapot')
title(['He:' num2str(He*1000) 'mV | Hi:' num2str(Hi*1000) 'mV'])


print -dpdf 'AnalAugust_He4k5Hi11k427' -r1000 
AnalAugust_He4k5Hi11k427={fixpointbranch,hopfbranch_psol_1,hopfbranch_psol_2}
save('AnalAugust_He4k5Hi11k427.mat','AnalAugust_He4k5Hi11k427')





% bif_RUM_default={fixpointbranch,hopfbranch_psol_1}
% save('bif_RUM_default.mat','bif_RUM_default')



% or just:
load('bif_RUM_default.mat')
fixpointbranch=bif_RUM_default{1}
hopfbranch_psol_1=bif_RUM_default{2}


% 
% load('bif_RUM_He34208Hi2200.mat')
% fixpointbranch=bif_2pop_He34208Hi2200{1}
% hopfbranch_psol_1=bif_2pop_He34208Hi2200{2}

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - Hopf branch - in He
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

% load('bif_RUM_default.mat')
% fixpointbranch=bif_RUM_default{1}
% hopfbranch_psol_1=bif_RUM_default{2}

ind_pext
ind_He=1;

hopfbif_branch_1=df_brnch(funcs,[ind_pext,ind_He],'hopf'); % use hopf point as first point of hopf branch:

pext_min= -80.0;
pext_max = 250.0;
pext_step = 1.0;


He_min = 1.5e-3;
He_max = 4.5e-3;
He_step = 0.05e-3;

ind
indexHopf       = 426;
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_1.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_He He_min]']';
hopfbif_branch_1.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_He He_max]']';
hopfbif_branch_1.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_He He_step]']';
hopfbif_branch_1.point=hopfpoint;
hopfbif_branch_1.method.continuation.plot=1;


hopfpoint.parameter(ind_He)=hopfpoint.parameter(ind_He)-0.001e-5; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_1.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,80)           % continue with plotting hopf branch:
hopfbif_branch_1=br_rvers(hopfbif_branch_1)                             % reverse Hopf branch
[hopfbif_branch_1,s,f,r]=br_contn(funcs,hopfbif_branch_1,600)           % continue in other direction
xlabel('pext');ylabel('He');



n=1
for i=10:2:length(hopfbif_branch_1.point)
    floq(n,1)=hopfbif_branch_1.point(i).parameter(ind_pext);
    floq(n,2)=hopfbif_branch_1.point(i).parameter(ind_He);
    hopfpoint=hopfbif_branch_1.point(i);
    PamRange(1)=hopfpoint.parameter(ind_pext) - 0.2*abs(hopfpoint.parameter(ind_pext));
    PamRange(3)=hopfpoint.parameter(ind_pext) + 0.2*abs(hopfpoint.parameter(ind_pext));
    PamRange(2)=(abs(PamRange(3))-abs(PamRange(2)))/30
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

figure(3)
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
plot(bif_hopf_He(1,:), bif_hopf_He(2,:))

print -dpdf 'RUM_2PamBIF_Hepext_StabilityHopfBifBranch' -r1000 



bif_hopf_He = zeros(2,length(hopfbif_branch_1.point));
for i=1:length(hopfbif_branch_1.point)
    bif_hopf_He(1,i)=hopfbif_branch_1.point(i).parameter(ind_pext);
    bif_hopf_He(2,i)=hopfbif_branch_1.point(i).parameter(ind_He);
end

figure(3); clf();
plot(bif_hopf_He(1,:), bif_hopf_He(2,:))


figure(16); clf();
plot(bif_hopf_He(1,:), bif_hopf_He(2,:))
hold on;
plot([-100 150],[3.25e-3 3.25e-3],'r.-')
xlabel('pext');ylabel('He');
title('Variation in He - RUM')
legend('hopfbifurcation branch')



hopfbif_branch_1=br_stabl(funcs,hopfbif_branch_1,0,0);

figure(7); clf();
for i=1:length(hopfbif_branch_1.point)
    a(i)=hopfbif_branch_1.point(i).omega
    %plot(hopfbif_branch_1.point(i).omega)
end
plot(a)
a=hopfbif_branch_1.point(1:20).omega
[xm,ym]=df_measr(1,hopfbif_branch_1); % plot stability versus point number:
ym.subfield='l0';
br_plot(hopfbif_branch_1,[],ym,'c');
ym.subfield='l1';
br_plot(hopfbif_branch_1,[],ym,'b');
xlabel('point number along branch');ylabel('\Re(\lambda)');


% plot point to parameter mapping
figure();

subplot(211)
plot(1:length(hopfbif_branch_1.point),bif_pam(:,1))
xlabel('point number')
ylabel('pext')
subplot(212)
plot(1:length(hopfbif_branch_1.point),bif_pam(:,2))
xlabel('point number')
ylabel('He')
title('point - parameter mapping')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - lower fold - in He
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_pext
ind_He=1;

lowerfold_branch_He=df_brnch(funcs,[ind_pext,ind_He],'fold'); 

pext_min= -80.0;
pext_max = 235.51;
pext_step = 0.5;


He_min = 1.0e-3;
He_max = 4.5e-3;
He_step = 0.015e-3;

ind
indexLowerFold       = 285;
lowerfoldpoint=p_tofold(funcs,fixpointbranch.point(indexLowerFold));

lowerfold_branch_He.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_He He_min]']';
lowerfold_branch_He.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_He He_max]']';
lowerfold_branch_He.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_He He_step]']';
lowerfold_branch_He.point=lowerfoldpoint;
lowerfold_branch_He.method.continuation.plot=0;

lowerfoldpoint.parameter(ind_He)=lowerfoldpoint.parameter(ind_He)+0.001e-3; % perturb hopf point
[lowerfoldpoint,success]=p_correc(funcs,lowerfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
lowerfold_branch_He.point(2)=lowerfoldpoint;                                 % use as second point of hopf branch:


figure(25); clf;
[lowerfold_branch_He,s,f,r]=br_contn(funcs,lowerfold_branch_He,850)           % continue with plotting hopf branch:
lowerfold_branch_He=br_rvers(lowerfold_branch_He)                             % reverse Hopf branch
[lowerfold_branch_He,s,f,r]=br_contn(funcs,lowerfold_branch_He,800)           % continue in other direction
xlabel('pext');ylabel('He');

figure(26); clf();
for i=1:length(lowerfold_branch_He.point)
    bif_lfb_He(i,1)=lowerfold_branch_He.point(i).parameter(ind_pext);
    bif_lfb_He(i,2)=lowerfold_branch_He.point(i).parameter(ind_He);
end
plot(bif_lfb_He(:,1), bif_lfb_He(:,2))
hold on;
plot(bif_hopf_He(:,1), bif_hopf_He(:,2),'g')
plot([-100 150],[3.25e-3 3.25e-3],'r.-')
xlabel('pext');ylabel('He');
legend('fold bifurcation branch',' hopf bifurcation branch')
title('Variation in He - RUM')

lowerfold_branch_He=br_stabl(funcs,lowerfold_branch_He,0,0);

figure(16); clf;
[xm,ym]=df_measr(1,lowerfold_branch_He); % plot stability versus point number:
ym.subfield='l0';
br_plot(lowerfold_branch_He,[],ym,'c');
ym.subfield='l1';
br_plot(lowerfold_branch_He,[],ym,'b');
xlabel('point number along branch');ylabel('\Re(\lambda)');


% plot point to parameter mapping
figure();
subplot(211)
plot(1:length(lowerfold_branch_He.point),bif_lfb_He(:,1))
xlabel('point number')
ylabel('pext')
subplot(212)
plot(1:length(lowerfold_branch_He.point),bif_lfb_He(:,2))
xlabel('point number')
ylabel('He')
title('point - parameter mapping')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - upper fold - in He
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_pext
ind_He=1;

upperfold_branch=df_brnch(funcs,[ind_pext,ind_He],'fold'); 

pext_min= -80.0;
pext_max = 235.51;
pext_step = 0.3;


He_min = 1.7e-3;
He_max = 4.5e-3;
He_step = 0.01e-3;

ind
indexUpperFold       = 514;
upperfoldpoint=p_tofold(funcs,fixpointbranch.point(indexUpperFold));

upperfold_branch.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_He He_min]']';
upperfold_branch.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_He He_max]']';
upperfold_branch.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_He He_step]']';
upperfold_branch.point=upperfoldpoint;
upperfold_branch.method.continuation.plot=0;

upperfoldpoint.parameter(ind_He)=upperfoldpoint.parameter(ind_He)+0.001e-3; % perturb hopf point
[upperfoldpoint,success]=p_correc(funcs,upperfoldpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
upperfold_branch.point(2)=upperfoldpoint;                                 % use as second point of hopf branch:


figure(); clf;
[upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,550)           % continue with plotting hopf branch:
upperfold_branch=br_rvers(upperfold_branch)                             % reverse Hopf branch
[upperfold_branch,s,f,r]=br_contn(funcs,upperfold_branch,1100)           % continue in other direction
xlabel('pext');ylabel('He');

figure(15); clf();
for i=1:length(upperfold_branch.point)
    bif_ufb_He(i,1)=upperfold_branch.point(i).parameter(ind_pext);
    bif_ufb_He(i,2)=upperfold_branch.point(i).parameter(ind_He);
end
plot(bif_ufb_He(:,1), bif_ufb_He(:,2))
hold on;
plot(bif_lfb_He(:,1), bif_lfb_He(:,2),'k')
plot(bif_hopf_He(:,1), bif_hopf_He(:,2),'g')
plot([-100 150],[3.25e-3 3.25e-3],'r.-')
xlabel('pext');ylabel('He');
legend('upper fold bifurcation branch','lower fold bifurcation branch',' hopf bifurcation branch')
title('Variation in He - RUM')
print -dpdf 'RUM_2PamBIF_He_pext' -r1000 

upperfold_branch=br_stabl(funcs,upperfold_branch,0,0);

figure(16); clf;
[xm,ym]=df_measr(1,upperfold_branch); % plot stability versus point number:
ym.subfield='l0';
br_plot(upperfold_branch,[],ym,'c');
ym.subfield='l1';
br_plot(upperfold_branch,[],ym,'b');
xlabel('point number along branch');ylabel('\Re(\lambda)');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - homoclinic bif - in He
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_pext
ind_He=1;

homclinbif_branch=df_brnch(funcs,[ind_pext,ind_He],'hcli'); 

pext_min= -80.0;
pext_max = 150;
pext_step = 0.5;


He_min = 1.7e-3;
He_max = 4.5e-3;
He_step = 0.01e-3;


indexhomclinbif       = 83; %get index as final point in the hopfbranch
homClinpoint=p_tohcli(funcs,hopfbranch_psol_1.point(indexhomclinbif));

homclinbif_branch.parameter.min_bound(1:2,:)=[[ind_pext pext_min]' [ind_He He_min]']';
homclinbif_branch.parameter.max_bound(1:2,:)=[[ind_pext pext_max]' [ind_He He_max]']';
homclinbif_branch.parameter.max_step(1:2,:)=[[ind_pext pext_step]' [ind_He He_step]']';
homclinbif_branch.point=homClinpoint;
homclinbif_branch.method.continuation.plot=1;

homClinpoint.parameter(ind_He)=homClinpoint.parameter(ind_He)+0.0001e-3; % perturb hopf point
%[homClinpoint,success]=p_correc(funcs,homClinpoint,ind_pext,[],method.point); % correct hopf point, recompute stability
homclinbif_branch.point(2)=homClinpoint;                                 % use as second point of hopf branch:


figure(9); clf;
[homclinbif_branch,s,f,r]=br_contn(funcs,homclinbif_branch,30)           % continue with plotting hopf branch:
homclinbif_branch=br_rvers(homclinbif_branch)                             % reverse Hopf branch
[homclinbif_branch,s,f,r]=br_contn(funcs,homclinbif_branch,10)           % continue in other direction
xlabel('pext');ylabel('He');

figure(15); clf();
for i=1:length(homclinbif_branch.point)
    bif_hcb_He(i,1)=homclinbif_branch.point(i).parameter(ind_pext);
    bif_hcb_He(i,2)=homclinbif_branch.point(i).parameter(ind_He);
end
plot(bif_hcb_He(:,1), bif_hcb_He(:,2),'r')
hold on;
plot(bif_ufb_He(:,1), bif_ufb_He(:,2),'-.b')
plot(bif_lfb_He(:,1), bif_lfb_He(:,2),'b')
plot(bif_hopf_He(:,1), bif_hopf_He(:,2),'g')

plot([-100 150],[3.25e-3 3.25e-3],'k.-')
xlabel('pext');ylabel('He');
legend('homoclinic bifurcation branch','upper fold bifurcation branch','lower fold bifurcation branch',' hopf bifurcation branch')
title('Variation in He - RUM')
print -dpdf 'RUM_2PamBIF_He_pext' -r1000 


% bif_RUM_He_pext={bif_hfb_He,bif_hcb_He,bif_ufb_He,bif_lfb_He,bif_hopf_He}
% save('bif_RUM_He_pext.mat','bif_RUM_He_pext')


homclinbif_branch=br_stabl(funcs,homclinbif_branch,0,0);

figure(16); clf;
[xm,ym]=df_measr(1,homclinbif_branch); % plot stability versus point number:
ym.subfield='l0';
br_plot(homclinbif_branch,[],ym,'c');
ym.subfield='l1';
br_plot(homclinbif_branch,[],ym,'b');
xlabel('point number along branch');ylabel('\Re(\lambda)');

% plot point to parameter mapping
for i=1:length(hopfbranch_psol_1.point)
    hopfsol(i,1)=hopfbranch_psol_1.point(i).parameter(ind_pext);
end
figure()
plot(1:length(hopfbranch_psol_1.point),hopfsol(:,1))
xlabel('point number')
ylabel('pext')
title('point - parameter mapping')













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - hopf-fold bif - in He
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

disp('Find and continue fold of periodic orbits in pext and He');

ind_pext
ind_He=1;



% load('bif_RUM_He34208Hi2200.mat')
% fixpointbranch=bif_2pop_He34208Hi2200{1}
% hopfbranch_psol_1=bif_2pop_He34208Hi2200{2}



indexhopffoldbif       = 34; % get parameter as minimal val (fold) from hopfbranch_psol_1 (generate hopf fold for He=3.4208mV)
%hopfbranch_psol_1.point(34).parameter(6)
%initialize branch

pext_max = 50;
pext_step = 0.05;


hopfbranch_psol_1.parameter.max_step=[ind_pext,pext_max]; % remove step size restriction
[pfuncs,pbranch,suc]=SetupPOfold(funcs,hopfbranch_psol_1,indexhopffoldbif,...
    'contpar',[ind_pext,ind_He],'dir',ind_pext,'step',pext_step,'print_residual_info',1);
if suc
    disp('POFold initialization finished');
else
    warning('POFold initialization failed');
end

figure(2);clf()
pbranch=br_contn(pfuncs,pbranch,100);
   

bif_hfb_He = zeros(2,length(pbranch.point));
for i=1:length(pbranch.point)
    bif_hfb_He(1,i)=pbranch.point(i).parameter(ind_pext);
    bif_hfb_He(2,i)=pbranch.point(i).parameter(ind_He);
end

figure(3);clf();
plot(bif_hfb_He(1,:),bif_hfb_He(2,:),'y*')


figure(16); clf();
plot(bif_hfb_He(1,:),bif_hfb_He(2,:),'c')
hold on;
plot(bif_ufb_He(:,1), bif_ufb_He(:,2),'-.b')
plot(bif_lfb_He(:,1), bif_lfb_He(:,2),'b')
plot(bif_hopf_He(:,1), bif_hopf_He(:,2),'g')
plot(bif_hcb_He(:,1), bif_hcb_He(:,2),'r')


plot([-100 150],[3.25e-3 3.25e-3],'k.-')
xlabel('pext');ylabel('He');
legend('hopf fold bifurcation branch','upper fold bifurcation branch','lower fold bifurcation branch',' hopf bifurcation branch','homoclinic bifurcation branch')
title('Variation in He - RUM')

print -dpdf 'RUM_2PamBIF_He_pext' -r1000 


% bif_RUM_He_pext={bif_hfb_He,bif_hcb_He,bif_ufb_He,bif_lfb_He,bif_hopf_He}
% save('bif_RUM_He_pext.mat','bif_RUM_He_pext')

load('bif_RUM_He_pext.mat')
bif_hfb_He=bif_RUM_He_pext{1}
bif_hcb_He=bif_RUM_He_pext{2}
bif_ufb_He=bif_RUM_He_pext{3}
bif_lfb_He=bif_RUM_He_pext{4}
bif_hopf_He=bif_RUM_He_pext{5}











% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 
% %
% %
% %   Second Part: He as primary bifurcation parameter
% %
% %
% %
% %   INCLUDED IN [...]_3popRUM_HeHi
% %
% %
% %
% %
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% %% Prepare work space: Second Part: He as primary bifurcation parameter
% 
% run init_ddebiftool;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Definition of user functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% He=3.25e-3;             % V
% Hi=22e-3;               % V
% taue=0.01;              % s
% taui=0.02;              % s
% kp1=0;                  % connectivity of external noise to excitatory cells
% kp2=0;                  % connectivity of external noise to inhibitory cells
% kp3=0;                  % connectivity of external noise to pyramidal cells
% n_ext_ii=0;             % external noise to inhibitory cells
% n_ext_ei=0;             % external noise to excitatory cells
% n_ext_py=0;             % external noise to pyramidal cells
% pext=0;                % signal arriving from other areas
% b=taue/taui;            % ration of tau_e to tau_i
% b1=1;                   % Switch for EI: if 1, EI is included
% b2=1;                   % switch for input: if 0, input is fed to Py
% b3=1;                   % switch for self connectivity: if 1, selfconn is off
% Nep=135;                % Conn gain from py to ei
% Npe=0.8*Nep;            % Conn gain from ei to py
% Npi=0.25*Nep;           % Conn gain from ii to py
% Nip=0.25*Nep;           % Conn gain from py to ii
% Npp=0.95*Nip;           % self connectivity ei: assumption Npp~Nip
% Nii=0.95*Npi;           % self connectivity ii: assumption Nii~Npi
% e0=2.5;                 % 1/s
% r=560;                  % 1/V
% v0=6e-3;                % V
% 
% par = zeros(24,1);
% 
% par(1) = He;
% par(2) = Hi;
% par(3) = taue;
% par(4) = taui;
% par(5) = b;
% par(6) = pext;
% par(7) = Nep;
% par(8) = Npe;
% par(9) = kp1;            % ext to EI
% par(10) = kp2;            % ext to II
% par(11) = kp3;            % ext to Py
% par(12) = n_ext_ei;       % ext to EI
% par(13) = n_ext_ii;       % ext to II
% par(14) = n_ext_py;       % ext to Py
% par(15) = b1;           % Switch for EI: if 1, EI is included
% par(16) = b2;           % switch for input: if 0, input is fed to Py
% par(17) = b3;           % switch for self connectivity: if 1, selfconn is off
% par(18) = e0;
% par(19) = r;
% par(20) = v0;
% par(21) = Npi;
% par(22) = Nip;
% par(23) = Npp;
% par(24) = Nii;
% par(25) = 0.0001; 
% 
% 
% % %generic mapped equation system:
% % dy(1) = x(6);
% % dy(2) = x(7);
% % dy(3) = x(8);
% % dy(4) = x(9);
% % dy(5) = x(10);
% % dy(6) = He*taue*r*       (Nep*      sig(x(2)-x(3)) + kp1*n_ext_ei + b2*pext)                                 - 2*x(6)   - x(1);
% % dy(7) = He*taue*r*       (b1*Npe*   sig(x(1)     ) + kp3*n_ext_py +(1-b2)*pext + (1-b3)*Npp*sig(x(2)-x(3)) ) - 2*x(7)   - x(2);
% % dy(8) = Hi*taui*b^2*r*Npi*            sig(x(4)-x(5))                                                           - 2*b*x(8) - b^2*x(3);
% % dy(9) = He*taue*r*(Nip*sig(x(2)-x(3)) + kp2*n_ext_ii )                                           - 2*x(9)   - x(4);
% % dy(10)= Hi*taui*b^2*r*Nii*(1-b3)*sig(x(4)-x(5))                                              - 2*b*x(10) - b^2*x(5);
% % 
% % MappedModelGeneric_sys_rhs=@(x,par)[...
% %     x(6,1);...
% %     x(7,1);...
% %     x(8,1);...
% %     x(9,1);...
% %     x(10,1);...
% %     He*taue*r*       (Nep*      sig(x(2,1)-x(3,1)) + kp1*n_ext_ei + b2*pext)                                 - 2*x(6,1)   - x(1,1);...
% %     He*taue*r*       (b1*Npe*   sig(x(1,1)     ) + kp3*n_ext_py +(1-b2)*pext + (1-b3)*Npp*sig(x(2,1)-x(3,1)) ) - 2*x(7,1)   - x(2,1);...
% %     Hi*taui*b^2*r*Npi*            sig(x(4,1)-x(5,1))                                                           - 2*b*x(8,1) - b^2*x(3,1);...
% %     He*taue*r*(Nip*sig(x(2,1)-x(3,1)) + kp2*n_ext_ii )                                           - 2*x(9,1)   - x(4,1);...
% %     Hi*taui*b^2*r*Nii*(1-b3)*sig(x(4,1)-x(5,1))                                              - 2*b*x(10,1) - b^2*x(5,1)]
% 
% sig = @(v,par) 2*par(18) ./ (1+exp(par(19)*par(20))*exp(-v));
% 
% MappedModelGeneric_sys_rhs=@(x,par)[...
%     x(6,1);...
%     x(7,1);...
%     x(8,1);...
%     x(9,1);...
%     x(10,1);...
%     par(1)*par(3)*par(19)*    (par(7)*      sig(x(2,1)-x(3,1),par) + par(9)*par(12) + par(16)*par(6))                                                - 2*x(6,1)   - x(1,1);...
%     par(1)*par(3)*par(19)*    (par(15)*par(8)* sig(x(1,1),par) + par(11)*par(14) + (1-par(16))*par(6) + (1-par(17))*par(23)*sig(x(2,1)-x(3,1),par) ) - 2*x(7,1)   - x(2,1);...
%     par(2)*par(4)*par(5)^2*   par(19)*par(21)*            sig(x(4,1)-x(5,1),par)                                                                     - 2*par(5)*x(8,1) - par(5)^2*x(3,1);...
%     par(1)*par(3)*par(19)*    (par(22)*sig(x(2,1)-x(3,1),par) + par(10)*par(13))                                                                     - 2*x(9,1)   - x(4,1);...
%     par(2)*par(4)*par(5)^2*par(19)*par(24)*(1-par(17))    *sig(x(4,1)-x(5,1),par)                                                                    - 2*par(5)*x(10,1) - par(5)^2*x(5,1)];
% 
% 
% 
% 
% % Delays
% neuron_tau=@()[25];
% 
% % Bifurcation parameter
% ind_He=1;
% 
% funcs=set_funcs(...
%     'sys_rhs',MappedModelGeneric_sys_rhs,...
%     'sys_tau',neuron_tau)
% 
% %% Initialise the first guessed fixed point
% stst.kind='stst';
% stst.parameter=par';
% stst.x=zeros(10,1);
% %stst.x=[6.1090; 9.2369; 5.8561; 1.5517; -0.0; 0; 0; 0; 0; 0];    % Stele der Hopfbif in default RUM
% %stst.x=[7.3758; 9.6540; 6.1888; 1.6163; -0.0; 0; 0; 0; 0; 0];    %Stelle des Fixpoint an pext=50 im oberen branch f?r default RUM
% %% Calculate Fixpoint curve 
% 
% flag_newhheur=1; % use the new steplength heuristic (Verheyden_2007)
% method=df_mthod(funcs,'stst',flag_newhheur);
% method.stability.minimal_real_part=-40
% [stst,success] = p_correc(funcs,stst,[],[],method.point)
% stst.x
% if ~success,
%     fprintf(1,'correction failed\n')
% end
% 
% % compute its stability:
% stst.stability = p_stabil(funcs,stst,method.stability)
% method.stability.minimal_real_part=-100; 
% stst.stability=p_stabil(funcs,stst,method.stability); % recompute stability:
% figure(1); clf;
% p_splot(stst); 
% 
% 
% 
% 
% % Initialize branch of trivial equilibria
% % get an empty branch with ind_pext as a free parameter:
% fixpointbranch=df_brnch(funcs,ind_He,'stst')
% 
% % set bounds for continuation parameter
% freepammin= 1.5e-3;
% freepammax = 9e-3;
% freepamstep = 0.1e-3;
% 
% fixpointbranch.parameter.min_bound(1,:)=[ind_He freepammin];
% fixpointbranch.parameter.max_bound(1,:)=[ind_He freepammax];
% fixpointbranch.parameter.max_step(1,:)=[ind_He freepamstep];
% % use stst as a first branch point:
% fixpointbranch.point=stst;
% 
% 
% %  Extend and continue branch of trivial equilibria
% stst.parameter(ind_He)=stst.parameter(ind_He)+0.001e-3;
% method=df_mthod(funcs,'stst')
% [stst,success]=p_correc(funcs,stst,[],[],method.point)
% 
% % use as a second branch point:
% fixpointbranch.point(2)=stst;
% fixpointbranch.method.continuation.plot=1;     %switch off plotting
% fixpointbranch.method.stability.minimal_real_part=-100;
% 
% 
% figure(2);clf();
% % continue in one direction:
% [fixpointbranch,s,f,r]=br_contn(funcs,fixpointbranch,3300)
% % turn the branch around:
% fixpointbranch=br_rvers(fixpointbranch);
% % continue in the other direction:
% [fixpointbranch,s,f,r]=br_contn(funcs,fixpointbranch,300)
% fixpointbranch=br_rvers(fixpointbranch);
% 
% 
% % compute stability of the branch
% % the first 0 is how many points to skip between stab calculations
% % the second 0 is to not recalculate any stability already present
% fixpointbranch=br_stabl(funcs,fixpointbranch,0,0);
% 
% 
% % obtain suitable scalar measures to plot stability along branch:
% [xm,ym]=df_measr(1,fixpointbranch)
% figure(3); clf;
% br_plot(fixpointbranch,xm,ym,'b');     % plot stability along branch:
% ym.subfield='l0';
% br_plot(fixpointbranch,xm,ym,'c');     % l0: rightmost roots
% plot([freepammin freepammax],[0 0],'-.');
% xlabel('free parameter');ylabel('\Re\lambda');
% 
% 
% %Plot Eigenvalues for specific point
% ind_point=22;
% figure(4);clf();
% p_splot(fixpointbranch.point(ind_point));
% title(['point: ' num2str(ind_point) ' | parameter:' num2str(fixpointbranch.point(ind_point).parameter(ind_He))])
% 
% 
% % plot stability versus point number:
% figure(4); clf;
% br_plot(fixpointbranch,[],ym,'b');
% br_plot(fixpointbranch,[],ym,'b.');
% plot([0 30],[0 0],'-.');
% xlabel('point number along branch');ylabel('\Re(\lambda)');
% 
% % Plot all Eigenvalues along branch and along free parameter
% 
% branch_summar=zeros(22,length(fixpointbranch.point));
% for i=1:length(fixpointbranch.point)
%     branch_summar(1,i)=fixpointbranch.point(i).parameter(ind_He);   % parameter value
%     branch_summar(2,i)=i;                                         % point number
%     branch_summar(3:12,i)=fixpointbranch.point(i).x;                     % fix point vals
%     branch_summar(13:22,i)=fixpointbranch.point(i).stability.l0;         % stability of eigen values
% end
% 
% %Detect zero crossings
% ind=[[],[]];
% 
% for i=1:10
%     [ind_temp,points_temp]=fun_detectzerocross(real(branch_summar(12+i,:)),branch_summar(2,:),branch_summar(1,:));
%     ind=[ind [ind_temp;points_temp]];
%     %points=[points points_temp];
% end
% 
% figure(5);clf();
% for i=1:10
%     subplot(2,10,i)
%     plot(branch_summar(1,:),real(branch_summar(12+i,:))); 
%     hold on;
%     plot([branch_summar(1,1) branch_summar(1,end)], [0 0],'k.-');
%     subplot(2,10,10+i)
%     xlabel('input firing rate')
%     ylabel('\Re(\lambda)')
%     plot(branch_summar(2,:),real(branch_summar(12+i,:)));
%     hold on;
%     plot([branch_summar(1,1) branch_summar(1,end)], [0 0],'k.-');
%     xlabel('point number')
%     ylabel('\Re(\lambda)')
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   Continuation of first Hopf bifurcation in a periodic solution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ind
% 
% freePamRange(1) = 2.6e-3;
% freePamRange(2) = 0.05e-3;
% freePamRange(3) = 5e-3;
% indexHopf       = 148;       %get from eigenvalues plot or print ind and points
% 
% [ hopfbranch_psol_1 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_He,freePamRange );
% 
% figure(); clf();
% [xm,ym]=df_measr(1,hopfbranch_psol_1);
% ym.subfield='mu';
% br_plot(hopfbranch_psol_1,[],ym,'b');
% xlabel('point')
% ylabel('Re \lambda')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   Continuation of second Hopf bifurcation in a periodic solution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ind
% 
% freePamRange(1) = 2.5e-3;
% freePamRange(2) = 0.01e-3;
% freePamRange(3) = 3.5e-3;
% indexHopf       = 170;       %get from eigenvalues plot or print ind and points
% 
% [ hopfbranch_psol_2 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_He,freePamRange );
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   Put Branches together
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %first hopf branch
% [ extval_Hopf_1 ] = fun_prepareBranchPlot( hopfbranch_psol_1 , ind_He);
% [ extval_Hopf_2 ] = fun_prepareBranchPlot( hopfbranch_psol_2 , ind_He);
% 
% %fixpointbranch
% fix_pot = zeros(2,length(fixpointbranch.point));
% for i=1:length(fixpointbranch.point)
%     fix_pot(1,i)=fixpointbranch.point(i).parameter(ind_He);
%     fix_pot(2,i)=fixpointbranch.point(i).x(2) - fixpointbranch.point(i).x(3);
% end
% 
% figure(12);clf();
% plot(fix_pot(1,:),fix_pot(2,:)./560,'b')
% hold on;
% plot(extval_Hopf_1(1,:),extval_Hopf_1(2,:)./560,'r');
% plot(extval_Hopf_1(1,:),extval_Hopf_1(3,:)./560,'r');
% 
% plot(extval_Hopf_2(1,:),extval_Hopf_2(2,:)./560,'g');
% plot(extval_Hopf_2(1,:),extval_Hopf_2(3,:)./560,'g');
% 
% xlabel('He')
% ylabel('pyrapot')
% title(['pext:' num2str(pext) 's-1 | Hi:' num2str(Hi*1000) 'mV'])
% print -dpdf 'BIF_RUMhe_pext0_Hi22' -r1000 
% 
% bif_RUMhe_pext0_Hi22={fixpointbranch,hopfbranch_psol_1}
% save('bif_RUMhe_pext0_Hi22.mat','bif_RUMhe_pext0_Hi22')
% 
% % bif_RUM_default={fixpointbranch,hopfbranch_psol_1}
% % save('bif_RUM_default.mat','bif_RUM_default')
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
% %%   2 parameter bifurcation (He,Hi) - Hopf branch - in Hi
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% ind_He=1;
% ind_Hi=2;
% 
% 
% % load('bif_RUMhe_pext0_Hi22.mat')
% % fixpointbranch=bif_RUMhe_pext0_Hi22{1}
% % hopfbranch_psol_1=bif_RUMhe_pext0_Hi22{2}
% 
% 
% hopfbif_branch_HeHi=df_brnch(funcs,[ind_He,ind_Hi],'hopf'); % use hopf point as first point of hopf branch:
% 
% He_min= 2.5e-3;
% He_max = 15.0e-3;
% He_step = 0.05e-3;
% 
% 
% Hi_min = 10e-3;
% Hi_max = 30e-3;
% Hi_step = 0.2e-3;
% 
% ind
% indexHopf       = 148;
% hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));
% 
% hopfbif_branch_HeHi.parameter.min_bound(1:2,:)=[[ind_He He_min]' [ind_Hi Hi_min]']';
% hopfbif_branch_HeHi.parameter.max_bound(1:2,:)=[[ind_He He_max]' [ind_Hi Hi_max]']';
% hopfbif_branch_HeHi.parameter.max_step(1:2,:)=[[ind_He He_step]' [ind_Hi Hi_step]']';
% hopfbif_branch_HeHi.point=hopfpoint;
% hopfbif_branch_HeHi.method.continuation.plot=1;
% 
% 
% hopfpoint.parameter(ind_Hi)=hopfpoint.parameter(ind_Hi)+0.001e-3; % perturb hopf point
% [hopfpoint,success]=p_correc(funcs,hopfpoint,ind_He,[],method.point); % correct hopf point, recompute stability
% hopfbif_branch_HeHi.point(2)=hopfpoint;                                 % use as second point of hopf branch:
% 
% 
% figure(6); clf;
% [hopfbif_branch_HeHi,s,f,r]=br_contn(funcs,hopfbif_branch_HeHi,300)           % continue with plotting hopf branch:
% hopfbif_branch_HeHi=br_rvers(hopfbif_branch_HeHi)                             % reverse Hopf branch
% [hopfbif_branch_HeHi,s,f,r]=br_contn(funcs,hopfbif_branch_HeHi,30)           % continue in other direction
% xlabel('He');ylabel('Hi');
% 
% 
% 
% bif_hopf_HeHi = zeros(2,length(hopfbif_branch_HeHi.point));
% for i=1:length(hopfbif_branch_HeHi.point)
%     bif_hopf_HeHi(1,i)=hopfbif_branch_HeHi.point(i).parameter(ind_He);
%     bif_hopf_HeHi(2,i)=hopfbif_branch_HeHi.point(i).parameter(ind_Hi);
% end
% 
% figure(14); clf();
% plot(bif_hopf_HeHi(1,:), bif_hopf_HeHi(2,:))
% hold on;
% plot([2e-3 16e-3],[22e-3 22e-3],'r.-')
% plot([3.25e-3 3.25e-3],[10e-3 24e-3],'r.-')
% xlabel('He');ylabel('Hi');
% legend('hopfbifurcation branch');
% title('He-Hi- bifurcation continuation')
% 
% 
% 
% hopfbif_branch_HeHi=br_stabl(funcs,hopfbif_branch_HeHi,0,0);
% 
% figure(7); clf;
% [xm,ym]=df_measr(1,hopfbif_branch_HeHi); % plot stability versus point number:
% ym.subfield='l0';
% br_plot(hopfbif_branch_HeHi,[],ym,'c');
% ym.subfield='l1';
% br_plot(hopfbif_branch_HeHi,[],ym,'b');
% xlabel('point number along branch');ylabel('\Re(\lambda)');
% 
% 
% % plot point to parameter mapping
% figure();
% 
% subplot(211)
% plot(1:length(hopfbif_branch_HeHi.point),bif_pam(:,1))
% xlabel('point number')
% ylabel('pext')
% subplot(212)
% plot(1:length(hopfbif_branch_HeHi.point),bif_pam(:,2))
% xlabel('point number')
% ylabel('He')
% title('point - parameter mapping')
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   2 parameter bifurcation - lower fold - in Hi
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% ind_He=1;
% ind_Hi=2;
% 
% 
% % load('bif_RUMhe_pext0_Hi22.mat')
% % fixpointbranch=bif_RUMhe_pext0_Hi22{1}
% % hopfbranch_psol_1=bif_RUMhe_pext0_Hi22{2}
% 
% 
% lowerfold_branch_HeHi=df_brnch(funcs,[ind_He,ind_Hi],'fold'); 
% He_min= 2e-3;
% He_max = 15.0e-3;
% He_step = 0.05e-3;
% 
% 
% Hi_min = 5e-3;
% Hi_max = 30e-3;
% Hi_step = 0.2e-3;
% 
% ind
% indexLowerFold_HeHi       = 101;
% lowerfoldpoint=p_tofold(funcs,fixpointbranch.point(indexLowerFold_HeHi));
% 
% lowerfold_branch_HeHi.parameter.min_bound(1:2,:)=[[ind_He He_min]' [ind_Hi Hi_min]']';
% lowerfold_branch_HeHi.parameter.max_bound(1:2,:)=[[ind_He He_max]' [ind_Hi Hi_max]']';
% lowerfold_branch_HeHi.parameter.max_step(1:2,:)=[[ind_He He_step]' [ind_Hi Hi_step]']';
% lowerfold_branch_HeHi.point=lowerfoldpoint;
% lowerfold_branch_HeHi.method.continuation.plot=1;
% 
% lowerfoldpoint.parameter(ind_Hi)=lowerfoldpoint.parameter(ind_Hi)+0.00001e-3; % perturb hopf point
% [lowerfoldpoint,success]=p_correc(funcs,lowerfoldpoint,ind_He,[],method.point); % correct hopf point, recompute stability
% lowerfold_branch_HeHi.point(2)=lowerfoldpoint;                                 % use as second point of hopf branch:
% 
% 
% figure(25); clf;
% [lowerfold_branch_HeHi,s,f,r]=br_contn(funcs,lowerfold_branch_HeHi,850)           % continue with plotting hopf branch:
% lowerfold_branch_HeHi=br_rvers(lowerfold_branch_HeHi)                             % reverse Hopf branch
% [lowerfold_branch_HeHi,s,f,r]=br_contn(funcs,lowerfold_branch_HeHi,800)           % continue in other direction
% xlabel('He');ylabel('Hi');
% 
% 
% 
% bif_lfb_HeHi = zeros(2,length(lowerfold_branch_HeHi.point));
% for i=1:length(lowerfold_branch_HeHi.point)
%     bif_lfb_HeHi(1,i)=lowerfold_branch_HeHi.point(i).parameter(ind_He);
%     bif_lfb_HeHi(2,i)=lowerfold_branch_HeHi.point(i).parameter(ind_Hi);
% end
% 
% figure(26); clf();
% plot(bif_lfb_HeHi(1,:), bif_lfb_HeHi(2,:))
% 
% 
% figure(14); clf();
% plot(bif_lfb_HeHi(:,1), bif_lfb_HeHi(:,2))
% hold on;
% plot(bif_hopf_HeHi(1,:), bif_hopf_HeHi(2,:))
% plot([2e-3 16e-3],[22e-3 22e-3],'k.-')
% plot([3.25e-3 3.25e-3],[10e-3 24e-3],'k.-')
% xlabel('He');ylabel('Hi');
% legend('lower fold bifurcation branch','hopfbifurcation branch');
% title(['He-Hi- bifurcation continuation, pext=' num2str(pext) 'pps' ])
% 
% 
% %print -dpdf 'RUM_2PamBIF_HeHi_pext' -r1000 
% %bif_RUM_HeHi={bif_lfb_HeHi,bif_hopf_HeHi}
% % save('bif_RUM_HeHi.mat','bif_RUM_HeHi')
% 
% 
% % load('bif_RUM_HeHi.mat');
% % bif_lfb_HeHi=bif_RUM_HeHi{1};
% % bif_hopf_HeHi=bif_RUM_HeHi{2};
% 
% 
% 
% 
% lowerfold_branch_HeHi=br_stabl(funcs,lowerfold_branch_HeHi,0,0);
% 
% figure(16); clf;
% [xm,ym]=df_measr(1,lowerfold_branch_HeHi); % plot stability versus point number:
% ym.subfield='l0';
% br_plot(lowerfold_branch_HeHi,[],ym,'c');
% ym.subfield='l1';
% br_plot(lowerfold_branch_HeHi,[],ym,'b');
% xlabel('point number along branch');ylabel('\Re(\lambda)');
% 
% 
% % plot point to parameter mapping
% figure();
% subplot(211)
% plot(1:length(lowerfold_branch_HeHi.point),bif_lfb_He(:,1))
% xlabel('point number')
% ylabel('pext')
% subplot(212)
% plot(1:length(lowerfold_branch_HeHi.point),bif_lfb_He(:,2))
% xlabel('point number')
% ylabel('He')
% title('point - parameter mapping')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%   2 parameter bifurcation - upper fold - in Hi
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first: identify relevant parameters and ranges through series of codim 1
% % bifurcation diagrams
% 
% 
% ind_He=1;
% ind_Hi=2;
% 
% 
% % load('bif_RUMhe_pext0_Hi22.mat')
% % fixpointbranch=bif_RUMhe_pext0_Hi22{1}
% % hopfbranch_psol_1=bif_RUMhe_pext0_Hi22{2}
% 
% 
% upperfold_branch_HeHi=df_brnch(funcs,[ind_He,ind_Hi],'fold'); 
% He_min= 1e-3;
% He_max = 15.0e-3;
% He_step = 0.05e-3;
% 
% 
% Hi_min = 5e-3;
% Hi_max = 30e-3;
% Hi_step = 0.2e-3;
% 
% ind
% indexUpperFold       = 150;
% upperfoldpoint=p_tofold(funcs,fixpointbranch.point(indexUpperFold));
% 
% upperfold_branch_HeHi.parameter.min_bound(1:2,:)=[[ind_He He_min]' [ind_Hi Hi_min]']';
% upperfold_branch_HeHi.parameter.max_bound(1:2,:)=[[ind_He He_max]' [ind_Hi Hi_max]']';
% upperfold_branch_HeHi.parameter.max_step(1:2,:)=[[ind_He He_step]' [ind_Hi Hi_step]']';
% upperfold_branch_HeHi.point=upperfoldpoint;
% upperfold_branch_HeHi.method.continuation.plot=1;
% 
% upperfoldpoint.parameter(ind_Hi)=upperfoldpoint.parameter(ind_Hi)+0.001e-3; % perturb hopf point
% [upperfoldpoint,success]=p_correc(funcs,upperfoldpoint,ind_He,[],method.point); % correct hopf point, recompute stability
% upperfold_branch_HeHi.point(2)=upperfoldpoint;                                 % use as second point of hopf branch:
% 
% 
% figure(); clf;
% [upperfold_branch_HeHi,s,f,r]=br_contn(funcs,upperfold_branch_HeHi,550)           % continue with plotting hopf branch:
% upperfold_branch_HeHi=br_rvers(upperfold_branch_HeHi)                             % reverse Hopf branch
% [upperfold_branch_HeHi,s,f,r]=br_contn(funcs,upperfold_branch_HeHi,1100)           % continue in other direction
% xlabel('He');ylabel('Hi');
% 
% 
% 
% bif_ufb_HeHi = zeros(2,length(upperfold_branch_HeHi.point));
% for i=1:length(upperfold_branch_HeHi.point)
%     bif_ufb_HeHi(i,1)=upperfold_branch_HeHi.point(i).parameter(ind_He);
%     bif_ufb_HeHi(i,2)=upperfold_branch_HeHi.point(i).parameter(ind_Hi);
% end
% 
% figure(15); clf();
% plot(bif_ufb_HeHi(:,1), bif_ufb_HeHi(:,2))
% 
% 
% figure(14); clf();
% plot(bif_ufb_HeHi(:,1), bif_ufb_HeHi(:,2),'-.b')
% hold on;
% plot(bif_lfb_HeHi(:,1), bif_lfb_HeHi(:,2),'b')
% plot(bif_hopf_HeHi(1,:), bif_hopf_HeHi(2,:),'g')
% plot([2e-3 16e-3],[22e-3 22e-3],'k.-')
% plot([3.25e-3 3.25e-3],[10e-3 24e-3],'k.-')
% xlabel('He');ylabel('Hi');
% legend('upper fold bifurcation branch','lower fold bifurcation branch','hopfbifurcation branch');
% title(['He-Hi- bifurcation continuation, pext=' num2str(pext) 'pps' ])
% 
% 
% %print -dpdf 'RUM_2PamBIF_HeHi_pext' -r1000 
% %bif_RUM_HeHi={bif_ufb_HeHi,bif_lfb_HeHi,bif_hopf_HeHi}
% % save('bif_RUM_HeHi.mat','bif_RUM_HeHi')
% 
% 
% % load('bif_RUM_HeHi.mat');
% % bif_ufb_HeHi=bif_RUM_HeHi{1};
% % bif_lfb_HeHi=bif_RUM_HeHi{2};
% % bif_hopf_HeHi=bif_RUM_HeHi{3};
% 
% 
% 
% upperfold_branch_HeHi=br_stabl(funcs,upperfold_branch_HeHi,0,0);
% 
% figure(16); clf;
% [xm,ym]=df_measr(1,upperfold_branch_HeHi); % plot stability versus point number:
% ym.subfield='l0';
% br_plot(upperfold_branch_HeHi,[],ym,'c');
% ym.subfield='l1';
% br_plot(upperfold_branch_HeHi,[],ym,'b');
% xlabel('point number along branch');ylabel('\Re(\lambda)');
% 





