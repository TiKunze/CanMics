

%%  Hi as primary bifurcation parameter, He as secondary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%
%
% ATTENTION: EQUATIONS WERE MODIFIED, SO THAT B1 tunes EIN and external input
%
%            b1 stays b1   
%            b2 gets  b1
%            Nep is multiplied by b1
%            b3 at Npp gets b1

run init_ddebiftool;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of user functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
He=3.25e-3;             % V
Hi=42e-3;               % V
taue=0.01;              % s
taui=0.02;              % s
kp1=0;                  % connectivity of external noise to excitatory cells
kp2=0;                  % connectivity of external noise to inhibitory cells
kp3=0;                  % connectivity of external noise to pyramidal cells
n_ext_ii=0;             % external noise to inhibitory cells
n_ext_ei=0;             % external noise to excitatory cells
n_ext_py=0;             % external noise to pyramidal cells
pext=100;                % signal arriving from other areas
b=taue/taui;            % ration of tau_e to tau_i
b1=0.0;                   % Switch for EI: if 1, EI is included
b2=1e10;                % not important in this system
b3=0.0;                   % switch for self connectivity of Nii: if 1, selfconn is off
Nep=135;                % Conn gain from py to ei
Npe=1e3*Nep;            % Conn gain from ei to py
Npi=0.25*Nep;           % Conn gain from ii to py
Nip=0.25*Nep;           % Conn gain from py to ii
Npp=113.4%0.7407*Nep;          % self connectivity ei: assumption Npp~Nip
Nii=33.25%100000.95*Npi;           % self connectivity ii: assumption Nii~Npi
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
% dy(4) = x(9);
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

ind_He=1
ind_Hi=2
%% Initialise the first guessed fixed point
stst.kind='stst';
stst.parameter=par';
stst.x=zeros(10,1);
%stst.x=[6.1090; 9.2369; 5.8561; 1.5517; -0.0; 0; 0; 0; 0; 0];    % Stele der Hopfbif in default RUM
%stst.x=[7.3758; 9.6540; 6.1888; 1.6163; -0.0; 0; 0; 0; 0; 0];    %Stelle des Fixpoint an pext=50 im oberen branch f?r default RUM
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
fixpointbranch=df_brnch(funcs,ind_Hi,'stst')

% set bounds for continuation parameter
freepammin= -5.0e-3;
freepammax = 1600e-3;
freepamstep = 7.1e-3;

fixpointbranch.parameter.min_bound(1,:)=[ind_Hi freepammin];
fixpointbranch.parameter.max_bound(1,:)=[ind_Hi freepammax];
fixpointbranch.parameter.max_step(1,:)=[ind_Hi freepamstep];
% use stst as a first branch point:
fixpointbranch.point=stst;


%  Extend and continue branch of trivial equilibria
stst.parameter(ind_Hi)=stst.parameter(ind_Hi)+0.001e-3;
method=df_mthod(funcs,'stst')
[stst,success]=p_correc(funcs,stst,[],[],method.point)

% use as a second branch point:
fixpointbranch.point(2)=stst;
fixpointbranch.method.continuation.plot=0;     %switch off plotting
fixpointbranch.method.stability.minimal_real_part=-100;


figure(1);clf();
% continue in one direction:
[fixpointbranch,s,f,r]=br_contn(funcs,fixpointbranch,3000)
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
ind_point=22;
figure(1);clf();
p_splot(fixpointbranch.point(ind_point));
title(['point: ' num2str(ind_point) ' | parameter:' num2str(fixpointbranch.point(ind_point).parameter(ind_Hi))])


% plot stability versus point number:
figure(1); clf;
br_plot(fixpointbranch,[],ym,'b');
br_plot(fixpointbranch,[],ym,'b.');
plot([0 30],[0 0],'-.');
xlabel('point number along branch');ylabel('\Re(\lambda)');

% Plot all Eigenvalues along branch and along free parameter

branch_summar=zeros(22,length(fixpointbranch.point));
for i=1:length(fixpointbranch.point)
    branch_summar(1,i)=fixpointbranch.point(i).parameter(ind_Hi);   % parameter value
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

freePamRange(1) = -10.5e-3;
freePamRange(2) = 0.3e-3;
freePamRange(3) = 26e-3;
indexHopf       = 124;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_1 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_Hi,freePamRange,100,1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Continuation of second Hopf bifurcation in a periodic solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind

freePamRange(1) = -7.5e-3;
freePamRange(2) = 0.1e-3;
freePamRange(3) = 10.5e-3;
indexHopf       = 210;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_2 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_Hi,freePamRange,80,1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Put Branches together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first hopf branch
[ extval_Hopf_1 ] = fun_prepareBranchPlot( hopfbranch_psol_1 , ind_Hi);
[ extval_Hopf_2 ] = fun_prepareBranchPlot( hopfbranch_psol_2 , ind_Hi);



figure(14);clf();
%plot(fix_pot(1,:),fix_pot(2,:)./560,'b')
fun_plotFPcurve(fixpointbranch,ind,ind_Hi,14)
hold on;
%fun_plotPersol(hopfbranch_1, ind_pext)

plot(extval_Hopf_1(1,:),extval_Hopf_1(2,:)./560,'r');
plot(extval_Hopf_1(1,:),extval_Hopf_1(3,:)./560,'r');

plot(extval_Hopf_2(1,:),extval_Hopf_2(2,:)./560,'g');
plot(extval_Hopf_2(1,:),extval_Hopf_2(3,:)./560,'g');

xlabel('Hi')
ylabel('pyrapot')
title(['He:' num2str(He*1000) 'mV | Hi:' num2str(Hi*1000) 'mV | Nii:' num2str(Nii) '| pext:' num2str(pext)])

print -dpdf 'BIF_RUMhi_pext100' -r1000 
bif_RUMhi_pext100={fixpointbranch}
save('bif_RUMhi_pext100.mat','bif_RUMhi_pext100')

,hopfbranch_psol_1
% bif_RUM_default={fixpointbranch,hopfbranch_psol_1}
% save('bif_RUM_default.mat','bif_RUM_default')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation (He,Hi) - Hopf branch - in Hi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_He=1;
ind_Hi=2;

% 
% load('bif_RUMhi_pext0.mat')
% fixpointbranch=bif_RUMhi_pext0{1}
% hopfbranch_psol_1=bif_RUMhi_pext0{2}


hopfbif_branch_HiHe=df_brnch(funcs,[ind_Hi,ind_He],'hopf'); % use hopf point as first point of hopf branch:


Hi_min = 1e-3;
Hi_max = 30e-3;
Hi_step = 0.2e-3;


He_min= 2.5e-3;
He_max = 15.0e-3;
He_step = 0.05e-3;




ind
indexHopf       = 110;
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_HiHe.parameter.min_bound(1:2,:)=[[ind_Hi Hi_min]' [ind_He He_min]']';
hopfbif_branch_HiHe.parameter.max_bound(1:2,:)=[[ind_Hi Hi_max]' [ind_He He_max]']';
hopfbif_branch_HiHe.parameter.max_step(1:2,:)=[[ind_Hi Hi_step]' [ind_He He_step]']';
hopfbif_branch_HiHe.point=hopfpoint;
hopfbif_branch_HiHe.method.continuation.plot=1;


hopfpoint.parameter(ind_He)=hopfpoint.parameter(ind_He)+0.001e-3; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_Hi,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_HiHe.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_HiHe,s,f,r]=br_contn(funcs,hopfbif_branch_HiHe,100)           % continue with plotting hopf branch:
hopfbif_branch_HiHe=br_rvers(hopfbif_branch_HiHe)                             % reverse Hopf branch
[hopfbif_branch_HiHe,s,f,r]=br_contn(funcs,hopfbif_branch_HiHe,100)           % continue in other direction
xlabel('Hi');ylabel('He');

%%

n=1
for i=10:5:length(hopfbif_branch_HiHe.point)
    i
    floq(n,1)=hopfbif_branch_HiHe.point(i).parameter(ind_Hi);
    floq(n,2)=hopfbif_branch_HiHe.point(i).parameter(ind_He);
    hopfpoint=hopfbif_branch_HiHe.point(i);
    PamRange(1)=hopfpoint.parameter(ind_Hi) - 0.2*abs(hopfpoint.parameter(ind_Hi));
    PamRange(3)=hopfpoint.parameter(ind_Hi) + 0.2*abs(hopfpoint.parameter(ind_Hi));
    PamRange(2)=(abs(PamRange(3))-abs(PamRange(2)))/30;
    temp=fun_determ_Floquet_HpfBifBranch( funcs,hopfpoint,fixpointbranch,indexHopf,ind_Hi,PamRange,20)
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



bif_hopf_HiHe = zeros(2,length(hopfbif_branch_HiHe.point));
for i=1:length(hopfbif_branch_HiHe.point)
    bif_hopf_HiHe(1,i)=hopfbif_branch_HiHe.point(i).parameter(ind_Hi);
    bif_hopf_HiHe(2,i)=hopfbif_branch_HiHe.point(i).parameter(ind_He);
end

figure(3)
clf()
plot(bif_hopf_HiHe(1,:), bif_hopf_HiHe(2,:))



figure(16); clf();
plot(bif_hopf_HiHe(1,:), bif_hopf_HiHe(2,:))
hold on;
plot([3e-3 6.5e-3],[3.25e-3 3.25e-3],'r.-')
%plot([22e-3 22e-3],[2.5e-3 7e-3],'r.-')
xlabel('Hi');ylabel('He');
legend('hopfbifurcation branch');
title('Hi-He- bifurcation continuation')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - lower fold - in He
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_He=1;
ind_Hi=2;


% load('bif_RUMhi_pext0.mat')
% fixpointbranch=bif_RUMhi_pext0{1}
% hopfbranch_psol_1=bif_RUMhi_pext0{2}


lowerfold_branch_HiHe=df_brnch(funcs,[ind_Hi,ind_He],'fold'); 
He_min= 1e-3;
He_max = 25.0e-3;
He_step = 0.15e-3;


Hi_min = -5e-3;
Hi_max = 30e-3;
Hi_step = 0.2e-3;

ind
indexLowerFold_HiHe       = 210;
lowerfoldpoint=p_tofold(funcs,fixpointbranch.point(indexLowerFold_HiHe));

lowerfold_branch_HiHe.parameter.min_bound(1:2,:)=[[ind_Hi Hi_min]' [ind_He He_min]']';
lowerfold_branch_HiHe.parameter.max_bound(1:2,:)=[[ind_Hi Hi_max]' [ind_He He_max]']';
lowerfold_branch_HiHe.parameter.max_step(1:2,:)=[[ind_Hi Hi_step]' [ind_He He_step]']';
lowerfold_branch_HiHe.point=lowerfoldpoint;
lowerfold_branch_HiHe.method.continuation.plot=1;

lowerfoldpoint.parameter(ind_He)=lowerfoldpoint.parameter(ind_He)-0.01e-3; % perturb hopf point
[lowerfoldpoint,success]=p_correc(funcs,lowerfoldpoint,ind_Hi,[],method.point); % correct hopf point, recompute stability
lowerfold_branch_HiHe.point(2)=lowerfoldpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[lowerfold_branch_HiHe,s,f,r]=br_contn(funcs,lowerfold_branch_HiHe,250)           % continue with plotting hopf branch:
lowerfold_branch_HiHe=br_rvers(lowerfold_branch_HiHe)                             % reverse Hopf branch
[lowerfold_branch_HiHe,s,f,r]=br_contn(funcs,lowerfold_branch_HiHe,250)           % continue in other direction
xlabel('Hi');ylabel('He');

%%

bif_lfb_HeHi = zeros(2,length(lowerfold_branch_HiHe.point));
for i=1:length(lowerfold_branch_HiHe.point)
    bif_lfb_HeHi(1,i)=lowerfold_branch_HiHe.point(i).parameter(ind_Hi);
    bif_lfb_HeHi(2,i)=lowerfold_branch_HiHe.point(i).parameter(ind_He);
end

figure(3); clf();
plot(bif_lfb_HeHi(1,:), bif_lfb_HeHi(2,:))





figure(16); clf();
plot(bif_lfb_HeHi(1,:), bif_lfb_HeHi(2,:),'b')
hold on;
plot(bif_hopf_HiHe(1,:), bif_hopf_HiHe(2,:))
plot([3e-3 6.5e-3],[3.25e-3 3.25e-3],'r.-')
xlabel('Hi');ylabel('He');
legend('lower fold bifurcation','hopfbifurcation branch');
title('Hi-He- bifurcation continuation')



%print -dpdf 'RUM_2PamBIF_HiHe' -r1000 
%bif_RUM_HiHe={bif_lfb_HeHi,bif_hopf_HiHe,floq}
% save('bif_RUM_HiHe.mat','bif_RUM_HiHe')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - upper fold - in He
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams


ind_He=1;
ind_Hi=2;



% load('bif_RUMhi_pext0.mat')
% fixpointbranch=bif_RUMhi_pext0{1}
% hopfbranch_psol_1=bif_RUMhi_pext0{2}


upperfold_branch_HiHe=df_brnch(funcs,[ind_Hi,ind_He],'fold'); 
He_min= -1e-3;
He_max = 15.0e-3;
He_step = 0.05e-3;


Hi_min = -5e-3;
Hi_max = 30e-3;
Hi_step = 0.12e-3;

ind
indexUpperFold       = 120;
upperfoldpoint=p_tofold(funcs,fixpointbranch.point(indexUpperFold));

upperfold_branch_HiHe.parameter.min_bound(1:2,:)=[[ind_Hi Hi_min]' [ind_He He_min]']';
upperfold_branch_HiHe.parameter.max_bound(1:2,:)=[[ind_Hi Hi_max]' [ind_He He_max]']';
upperfold_branch_HiHe.parameter.max_step(1:2,:)=[[ind_Hi Hi_step]' [ind_He He_step]']';
upperfold_branch_HiHe.point=upperfoldpoint;
upperfold_branch_HiHe.method.continuation.plot=1;

upperfoldpoint.parameter(ind_He)=upperfoldpoint.parameter(ind_He)-0.001e-3; % perturb hopf point
[upperfoldpoint,success]=p_correc(funcs,upperfoldpoint,ind_Hi,[],method.point); % correct hopf point, recompute stability
upperfold_branch_HiHe.point(2)=upperfoldpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[upperfold_branch_HiHe,s,f,r]=br_contn(funcs,upperfold_branch_HiHe,150)           % continue with plotting hopf branch:
upperfold_branch_HiHe=br_rvers(upperfold_branch_HiHe)                             % reverse Hopf branch
[upperfold_branch_HiHe,s,f,r]=br_contn(funcs,upperfold_branch_HiHe,400)           % continue in other direction
xlabel('Hi');ylabel('He');

%%

bif_ufb_HiHe = zeros(2,length(upperfold_branch_HiHe.point));
for i=1:length(upperfold_branch_HiHe.point)
    bif_ufb_HiHe(1,i)=upperfold_branch_HiHe.point(i).parameter(ind_Hi);
    bif_ufb_HiHe(2,i)=upperfold_branch_HiHe.point(i).parameter(ind_He);
end

figure(3); clf();
plot(bif_ufb_HiHe(1,:), bif_ufb_HiHe(2,:))



figure(16); clf();
plot(bif_ufb_HiHe(1,:), bif_ufb_HiHe(2,:),'-.b')
hold on;
plot(bif_lfb_HeHi(1,:), bif_lfb_HeHi(2,:),'b')
plot(bif_hopf_HiHe(1,:), bif_hopf_HiHe(2,:))
plot([3e-3 6.5e-3],[3.25e-3 3.25e-3],'r.-')
xlabel('Hi');ylabel('He');
legend('upper fold bifurcation','lower fold bifurcation','hopfbifurcation branch');
title('Hi-He- bifurcation continuation')



%print -dpdf 'RUM_2PamBIF_HiHe' -r1000 
%bif_RUM_HiHe={bif_ufb_HiHe,bif_lfb_HeHi,bif_hopf_HiHe,floq}
% save('bif_RUM_HiHe.mat','bif_RUM_HiHe')








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation (He,Hi) - Hopf branch - in Hi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_He=1;
ind_Hi=2;

% 
% load('bif_RUMhi_pext0.mat')
% fixpointbranch=bif_RUMhi_pext0{1}
% hopfbranch_psol_1=bif_RUMhi_pext0{2}


hopfbif_branch_HiHe=df_brnch(funcs,[ind_Hi,ind_He],'hopf'); % use hopf point as first point of hopf branch:


Hi_min = 1e-3;
Hi_max = 30e-3;
Hi_step = 0.2e-3;


He_min= 2.5e-3;
He_max = 15.0e-3;
He_step = 0.05e-3;




ind
indexHopf       = 110;
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_HiHe.parameter.min_bound(1:2,:)=[[ind_Hi Hi_min]' [ind_He He_min]']';
hopfbif_branch_HiHe.parameter.max_bound(1:2,:)=[[ind_Hi Hi_max]' [ind_He He_max]']';
hopfbif_branch_HiHe.parameter.max_step(1:2,:)=[[ind_Hi Hi_step]' [ind_He He_step]']';
hopfbif_branch_HiHe.point=hopfpoint;
hopfbif_branch_HiHe.method.continuation.plot=1;


hopfpoint.parameter(ind_He)=hopfpoint.parameter(ind_He)+0.001e-3; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_Hi,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_HiHe.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_HiHe,s,f,r]=br_contn(funcs,hopfbif_branch_HiHe,100)           % continue with plotting hopf branch:
hopfbif_branch_HiHe=br_rvers(hopfbif_branch_HiHe)                             % reverse Hopf branch
[hopfbif_branch_HiHe,s,f,r]=br_contn(funcs,hopfbif_branch_HiHe,100)           % continue in other direction
xlabel('Hi');ylabel('He');

%%

n=1
for i=10:5:length(hopfbif_branch_HiHe.point)
    i
    floq(n,1)=hopfbif_branch_HiHe.point(i).parameter(ind_Hi);
    floq(n,2)=hopfbif_branch_HiHe.point(i).parameter(ind_He);
    hopfpoint=hopfbif_branch_HiHe.point(i);
    PamRange(1)=hopfpoint.parameter(ind_Hi) - 0.2*abs(hopfpoint.parameter(ind_Hi));
    PamRange(3)=hopfpoint.parameter(ind_Hi) + 0.2*abs(hopfpoint.parameter(ind_Hi));
    PamRange(2)=(abs(PamRange(3))-abs(PamRange(2)))/30;
    temp=fun_determ_Floquet_HpfBifBranch( funcs,hopfpoint,fixpointbranch,indexHopf,ind_Hi,PamRange,20)
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



bif_hopf_HiHe = zeros(2,length(hopfbif_branch_HiHe.point));
for i=1:length(hopfbif_branch_HiHe.point)
    bif_hopf_HiHe(1,i)=hopfbif_branch_HiHe.point(i).parameter(ind_Hi);
    bif_hopf_HiHe(2,i)=hopfbif_branch_HiHe.point(i).parameter(ind_He);
end

figure(3)
clf()
plot(bif_hopf_HiHe(1,:), bif_hopf_HiHe(2,:))



figure(16); clf();
plot(bif_hopf_HiHe(1,:), bif_hopf_HiHe(2,:))
hold on;
plot([3e-3 6.5e-3],[3.25e-3 3.25e-3],'r.-')
%plot([22e-3 22e-3],[2.5e-3 7e-3],'r.-')
xlabel('Hi');ylabel('He');
legend('hopfbifurcation branch');
title('Hi-He- bifurcation continuation')




