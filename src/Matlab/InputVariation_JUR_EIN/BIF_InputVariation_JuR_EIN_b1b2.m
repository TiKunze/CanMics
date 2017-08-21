
%{
This model is a specific version of the canonical processing unit and is 
desgined to easily compare the behavior of the JuR (input to Py) and the EIN (input to EIN) model

The architectural parameters p_ext_Py (input to Py) and p_ext_EIN (input to EIN) tune the input distribution.
the architectural parameter b1 tunes the inclusion of the EIN (b1=0: no EIN)


%}


%%  He as primary bifurcation parameter, Hi as secondary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
b=taue/taui;            % ration of tau_e to tau_i

pext=300;               % maximum signal arriving from other areas (scaled by r1/2)


p_ext_Py=150.0;                   % input to Py 0->no input to Py
p_ext_EIN=0.0;                   % input to EIN 0-> no input to EIN
b1=1.0;                   % Switch for EIN: b1=0: no EIN, but direct feedback (over Npp)
b2=0.0;                   % controls disinhibition: 1->disinhbition (Nii) || 0->no disinhibition

Nep=135;                % Conn gain from py to ei
Npe=0.8*Nep;            % Conn gain from ei to py
Npi=0.25*Nep;           % Conn gain from ii to py
Nip=0.25*Nep;           % Conn gain from py to ii
Npp=113.4;              % self connectivity ei: from mapping
Nii=33.25;              % self connectivity ii:
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
par(15) = b1;           
par(16) = b2;           
par(17) = p_ext_Py;
par(18) = p_ext_EIN;
par(19) = e0;
par(20) = r;
par(21) = v0;
par(22) = Npi;
par(23) = Nip;
par(24) = Npp;
par(25) = Nii;
par(26) = 0.0001; 


% %generic mapped equation system:
% dy(1) = x(6)
% dy(2) = x(7)
% dy(3) = x(8)
% dy(4) = x(9)
% dy(5) = x(10)        
% dy(6) = self.H_e*self.tau_e*self.r *          (self.c_ep*self.sigm(x(2)-x(3)) + self.r2*(self.c_In_ex*impact_ex + noi_ex) )                 - 2*x(6) - x(1)
% dy(7) = self.H_e*self.tau_e*self.r *          (self.b1*self.c_pe*self.sigm(x(1)) + self.r1*(self.c_In_ex*impact_ex + noi_ex) + (1-self.b1)*self.c_pp*self.sigm(x(2)-x(3)) ) - 2*x(7) - x(2)  
% dy(8) = self.H_i*self.tau_i*self.eta**2*self.r * (self.c_pi * self.sigm(x(4)-x(5)))                                                                   - 2*self.eta*x(8) - self.eta**2*x(3)
% dy(9) = self.H_e*self.tau_e*self.r *          (self.c_In_ii*impact_i + noi_i + self.c_ip*self.sigm(x(2)-x(3)))                                             - 2*x(9) - x(4)  
% dy(10) = self.H_i*self.tau_i*self.eta**2*self.r * ((self.b2)*self.c_ii*self.sigm(x(4)-x(5)))                                                         - 2*self.eta*x(10) - self.eta**2*x(5)
% 
% dy(1) = x(6)
% dy(2) = x(7)
% dy(3) = x(8)
% dy(4) = x(9)
% dy(5) = x(10)        
% dy(6) = He*taue*r *          (Nep*sig(x(2)-x(3)) + p_ext_EIN )                 - 2*x(6) - x(1)
% dy(7) = He*taue*r *          (b1*par(8)*sig(x(1)) + p_ext_Py + (1-b1)*Npp*sig(x(2)-x(3)) ) - 2*x(7) - x(2)  
% dy(8) = Hi*taui*b^2*r * (Npi * sig(x(4)-x(5)))                                                                   - 2*b*x(8) - b^2*x(3)
% dy(9) = He*taue*r *          (Nip*sig(x(2)-x(3)))                                             - 2*x(9) - x(4)  
% dy(10) = Hi*taui*b^2*r * b2*Nii*sig(x(4)-x(5))                                                       - 2*b*x(10) - b^2*x(5)
% 

sig = @(v,par) 2*par(19) ./ (1+exp(par(20)*par(21))*exp(-v));
MappedModelGeneric_sys_rhs=@(x,par)[...
    x(6,1);...
    x(7,1);...
    x(8,1);...
    x(9,1);...
    x(10,1);...        
    par(1)*par(3)*par(20) *          (par(7)*sig(x(2,1)-x(3,1),par) + par(18) )                 - 2*x(6,1) - x(1,1);...
    par(1)*par(3)*par(20) *          (par(15)*Npe*sig(x(1,1),par) + par(17) + (1-par(15))*par(24)*sig(x(2,1)-x(3,1),par) ) - 2*x(7,1) - x(2,1)  ;...
    par(2)*par(4)*par(5)^2*par(20) * (par(22) * sig(x(4,1)-x(5,1),par))                                                                   - 2*par(5)*x(8,1) - par(5)^2*x(3,1);...
    par(1)*par(3)*par(20) *          (par(23)*sig(x(2,1)-x(3,1),par))                                             - 2*x(9,1) - x(4,1)  ;...
    par(2)*par(4)*par(5)^2*par(20) * par(16)*par(25)*sig(x(4,1)-x(5,1),par)                                                       - 2*par(5)*x(10,1) - par(5)^2*x(5,1)];
 


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




% Delays
neuron_tau=@()[25];

% Bifurcation parameter
ind_pPy=17;
ind_pEIN=18;

funcs=set_funcs(...
    'sys_rhs',MappedModelGeneric_sys_rhs,...
    'sys_tau',neuron_tau)

% Initialise the first guessed fixed point
stst.kind='stst';
stst.parameter=par';
stst.x=zeros(10,1);

%% Calculate Fixpoint curve 

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
fixpointbranch=df_brnch(funcs,ind_pPy,'stst')

% set bounds for continuation parameter
freepammin= -100;
freepammax = 350;
freepamstep = 0.25; %1.5

fixpointbranch.parameter.min_bound(1,:)=[ind_pPy freepammin];
fixpointbranch.parameter.max_bound(1,:)=[ind_pPy freepammax];
fixpointbranch.parameter.max_step(1,:)=[ind_pPy freepamstep];
% use stst as a first branch point:
fixpointbranch.point=stst;


%  Extend and continue branch of trivial equilibria
stst.parameter(ind_pPy)=stst.parameter(ind_pPy)+0.001e-3;
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
[fixpointbranch,s,f,r]=br_contn(funcs,fixpointbranch,300)
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
title(['point: ' num2str(ind_point) ' | parameter:' num2str(fixpointbranch.point(ind_point).parameter(ind_pPy))])


% plot stability versus point number:
figure(1); clf;
br_plot(fixpointbranch,[],ym,'b');
br_plot(fixpointbranch,[],ym,'b.');
plot([0 30],[0 0],'-.');
xlabel('point number along branch');ylabel('\Re(\lambda)');

% Plot all Eigenvalues along branch and along free parameter

branch_summar=zeros(22,length(fixpointbranch.point));
for i=1:length(fixpointbranch.point)
    branch_summar(1,i)=fixpointbranch.point(i).parameter(ind_pPy);   % parameter value
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

freePamRange(1) = -15;
freePamRange(2) = 1.5;
freePamRange(3) = 150;
indexHopf       = 1572%452;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_1 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pPy,freePamRange,250,1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Continuation of second Hopf bifurcation in a periodic solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind

freePamRange(1) = 80;
freePamRange(2) = 2.5;
freePamRange(3) = 318;
indexHopf       = 1980%520;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_2 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pPy,freePamRange,150,1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Put Branches together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first hopf branch
[ extval_Hopf_1 ] = fun_prepareBranchPlot( hopfbranch_psol_1 , ind_pPy);
[ extval_Hopf_2 ] = fun_prepareBranchPlot( hopfbranch_psol_2 , ind_pPy);



figure(12);clf();
%plot(fix_pot(1,:),fix_pot(2,:)./560,'b')
a=fun_plotFPcurve(fixpointbranch,ind,ind_pPy,12);
hold on;


plot(extval_Hopf_1(1,:),extval_Hopf_1(2,:)./560,'r');
plot(extval_Hopf_1(1,:),extval_Hopf_1(3,:)./560,'r');

plot(extval_Hopf_2(1,:),extval_Hopf_2(2,:)./560,'g');
plot(extval_Hopf_2(1,:),extval_Hopf_2(3,:)./560,'g');

xlabel('pext Py')
ylabel('pyrapot')
title(['pext EIN:' num2str(p_ext_EIN) 's-1 '])



print -dpdf 'BIF_InVar_EIN_0_JuR' -r1000 

bif_InVar_EIN_0_JuR={fixpointbranch,hopfbranch_psol_1,hopfbranch_psol_2}
save('bif_InVar_EIN_0_JuR.mat','bif_InVar_EIN_0_JuR')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation (p_Py,p_EIN) - subcritical Hopf branch - in pEIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_pEIN=18;
ind_pPy=17;


% load('BIF_InVar_EIN_0_JuR.mat')
% fixpointbranch=BIF_InVar_EIN_0_JuR{1}
% hopfbranch_psol_1=BIF_InVar_EIN_0_JuR{2}
% hopfbranch_psol_2=BIF_InVar_EIN_0_JuR{3}

hopfbif_branch_pPyEIN_1=df_brnch(funcs,[ind_pPy,ind_pEIN],'hopf'); % use hopf point as first point of hopf branch:

pPy_min= -100;
pPy_max = 350;
pPy_step = 2.5;


pEIN_min = -100;
pEIN_max = 350;
pEIN_step = 2.5;

ind
indexHopf       = 1572 %452;  520, 671
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_pPyEIN_1.parameter.min_bound(1:2,:)=[[ind_pPy pPy_min]' [ind_pEIN pEIN_min]']';
hopfbif_branch_pPyEIN_1.parameter.max_bound(1:2,:)=[[ind_pPy pPy_max]' [ind_pEIN pEIN_max]']';
hopfbif_branch_pPyEIN_1.parameter.max_step(1:2,:)=[[ind_pPy pPy_step]' [ind_pEIN pEIN_step]']';
hopfbif_branch_pPyEIN_1.point=hopfpoint;
hopfbif_branch_pPyEIN_1.method.continuation.plot=1;


hopfpoint.parameter(ind_pEIN)=hopfpoint.parameter(ind_pEIN)+0.001e-3; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pPy,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_pPyEIN_1.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_pPyEIN_1,s,f,r]=br_contn(funcs,hopfbif_branch_pPyEIN_1,170)           % continue with plotting hopf branch:
hopfbif_branch_pPyEIN_1=br_rvers(hopfbif_branch_pPyEIN_1)                             % reverse Hopf branch
[hopfbif_branch_pPyEIN_1,s,f,r]=br_contn(funcs,hopfbif_branch_pPyEIN_1,100)           % continue in other direction
xlabel('pPy');ylabel('pEIN');

%

bif_hopf_pPypEIN_1 = zeros(2,length(hopfbif_branch_pPyEIN_1.point));
for i=1:length(hopfbif_branch_pPyEIN_1.point)
    bif_hopf_pPypEIN_1(1,i)=hopfbif_branch_pPyEIN_1.point(i).parameter(ind_pPy);
    bif_hopf_pPypEIN_1(2,i)=hopfbif_branch_pPyEIN_1.point(i).parameter(ind_pEIN);
end

figure(3); clf();
plot(bif_hopf_pPypEIN_1(2,:), bif_hopf_pPypEIN_1(1,:))
hold on;

% in the following routine, (a short part of) the periodic solution for a certain point along the
% Hopf-bifurcation branch is calculated. A function fun_determ_Floquet_HpfBifBranch
% is called which provides Floquet-Multipliers. those are used to determine whether the periodic
% solution is stable or not. However, a stable point can rapidly change in
% to a unstable by means of a hopf-fold bifurcation
n=1
for i=1:4:length(hopfbif_branch_pPyEIN_1.point)
    i
    floq_hopf1(n,1)=hopfbif_branch_pPyEIN_1.point(i).parameter(ind_pEIN);   % extract current pEIN value
    floq_hopf1(n,2)=hopfbif_branch_pPyEIN_1.point(i).parameter(ind_pPy);   % extract current pPy value
    hopfpoint=hopfbif_branch_pPyEIN_1.point(i);                     % extract current data
    PamRange(1)=hopfpoint.parameter(ind_pEIN) - 0.2*abs(hopfpoint.parameter(ind_pEIN)); % estimate lower border of periodic solution 
    PamRange(3)=hopfpoint.parameter(ind_pEIN) + 0.2*abs(hopfpoint.parameter(ind_pEIN)); % estimate upper border of periodic solution
    PamRange(2)=(abs(PamRange(3))-abs(PamRange(2)))/30;
    temp=fun_determ_Floquet_HpfBifBranch( funcs,hopfpoint,fixpointbranch,indexHopf,ind_pEIN,PamRange,20)
    if length(temp)==2
        floq_hopf1(n,3:4)=temp;
    elseif length(temp)==4
        floq_hopf1(n,3:6)=temp;
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
for i=1:length(floq_hopf1(:,1))
    if max(floq_hopf1(i,3:4))>1.005
        plot(floq_hopf1(i,1),floq_hopf1(i,2),'.r','MarkerSize',20)
    elseif max(floq_hopf1(i,3:4))<1.001    
        plot(floq_hopf1(i,1),floq_hopf1(i,2),'.g','MarkerSize',20)
    else
        plot(floq_hopf1(i,1),floq_hopf1(i,2),'.m','MarkerSize',20)
    end
end





figure(16); clf();
plot(bif_hopf_pPypEIN_1(2,:), bif_hopf_pPypEIN_1(1,:))
hold on;
xlabel('pEIN');ylabel('pPy');
legend('hopfbifurcation branch');
title('pPy-pEIN- bifurcation continuation')


% 
% print -dpdf 'InVar_pPypEIN' -r1000 
% bif_InVar_pPypEIN={bif_hopf_pPypEIN,floq_hopf1}
% save('bif_InVar_pPypEIN.mat','bif_InVar_pPypEIN')



 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation (p_Py,p_EIN) - supercritical Hopf branch 2 - in pEIN
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schl?gt momentan fehl, aber ist vermutlich eh der gleiche Branch wie
% supercritical hpf 1

% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_pEIN=18;
ind_pPy=17;


% load('BIF_InVar_EIN_0_JuR.mat')
% fixpointbranch=BIF_InVar_EIN_0_JuR{1}
% hopfbranch_psol_1=BIF_InVar_EIN_0_JuR{2}
% hopfbranch_psol_2=BIF_InVar_EIN_0_JuR{3}

hopfbif_branch_pPyEIN_2=df_brnch(funcs,[ind_pPy,ind_pEIN],'hopf'); % use hopf point as first point of hopf branch:

pPy_min= -100;
pPy_max = 350;
pPy_step = 2.5;


pEIN_min = -100;
pEIN_max = 350;
pEIN_step = 2.5;

ind
indexHopf       = 1980%520 %, 671
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_pPyEIN_2.parameter.min_bound(1:2,:)=[[ind_pPy pPy_min]' [ind_pEIN pEIN_min]']';
hopfbif_branch_pPyEIN_2.parameter.max_bound(1:2,:)=[[ind_pPy pPy_max]' [ind_pEIN pEIN_max]']';
hopfbif_branch_pPyEIN_2.parameter.max_step(1:2,:)=[[ind_pPy pPy_step]' [ind_pEIN pEIN_step]']';
hopfbif_branch_pPyEIN_2.point=hopfpoint;
hopfbif_branch_pPyEIN_2.method.continuation.plot=1;


hopfpoint.parameter(ind_pEIN)=hopfpoint.parameter(ind_pEIN)+0.01; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pPy,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_pPyEIN_2.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_pPyEIN_2,s,f,r]=br_contn(funcs,hopfbif_branch_pPyEIN_2,250)           % continue with plotting hopf branch:
hopfbif_branch_pPyEIN_2=br_rvers(hopfbif_branch_pPyEIN_2)                             % reverse Hopf branch
[hopfbif_branch_pPyEIN_2,s,f,r]=br_contn(funcs,hopfbif_branch_pPyEIN_2,250)           % continue in other direction
xlabel('pPy');ylabel('pEIN');

%

bif_hopf_pPypEIN_2 = zeros(2,length(hopfbif_branch_pPyEIN_2.point));
for i=1:length(hopfbif_branch_pPyEIN_2.point)
    bif_hopf_pPypEIN_2(1,i)=hopfbif_branch_pPyEIN_2.point(i).parameter(ind_pPy);
    bif_hopf_pPypEIN_2(2,i)=hopfbif_branch_pPyEIN_2.point(i).parameter(ind_pEIN);
end

figure(3); clf();
plot(bif_hopf_pPypEIN_2(2,:), bif_hopf_pPypEIN_2(1,:))
hold on;

% in the following routine, (a short part of) the periodic solution for a certain point along the
% Hopf-bifurcation branch is calculated. A function fun_determ_Floquet_HpfBifBranch
% is called which provides Floquet-Multipliers. those are used to determine whether the periodic
% solution is stable or not. However, a stable point can rapidly change in
% to a unstable by means of a hopf-fold bifurcation
n=1
for i=1:4:length(hopfbif_branch_pPyEIN_2.point)
    i
    floq_hopf2(n,1)=hopfbif_branch_pPyEIN_2.point(i).parameter(ind_pEIN);   % extract current pEIN value
    floq_hopf2(n,2)=hopfbif_branch_pPyEIN_2.point(i).parameter(ind_pPy);   % extract current pPy value
    hopfpoint=hopfbif_branch_pPyEIN_2.point(i);                     % extract current data
    PamRange(1)=hopfpoint.parameter(ind_pEIN) - 0.2*abs(hopfpoint.parameter(ind_pEIN)); % estimate lower border of periodic solution 
    PamRange(3)=hopfpoint.parameter(ind_pEIN) + 0.2*abs(hopfpoint.parameter(ind_pEIN)); % estimate upper border of periodic solution
    PamRange(2)=(abs(PamRange(3))-abs(PamRange(2)))/30;
    temp=fun_determ_Floquet_HpfBifBranch( funcs,hopfpoint,fixpointbranch,indexHopf,ind_pEIN,PamRange,20)
    if length(temp)==2
        floq_hopf2(n,3:4)=temp;
    elseif length(temp)==4
        floq_hopf2(n,3:6)=temp;
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
for i=1:length(floq_hopf2(:,1))
    if max(floq_hopf2(i,3:4))>1.005
        plot(floq_hopf2(i,1),floq_hopf2(i,2),'.r','MarkerSize',20)
    elseif max(floq_hopf2(i,3:4))<1.001    
        plot(floq_hopf2(i,1),floq_hopf2(i,2),'.g','MarkerSize',20)
    else
        plot(floq_hopf2(i,1),floq_hopf2(i,2),'.m','MarkerSize',20)
    end
end





figure(16); clf();
plot(hopfbif_branch_pPyEIN_2(2,:), hopfbif_branch_pPyEIN_2(1,:))
plot(hopfbif_branch_pPyEIN_1(2,:), hopfbif_branch_pPyEIN_1(1,:))
hold on;
xlabel('pEIN');ylabel('pPy');
legend('hopfbifurcation branch');
title('pPy-pEIN- bifurcation continuation')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation (p_Py,p_EIN) - supercritical Hopf branch 3 - in pEIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_pEIN=18;
ind_pPy=17;


% load('BIF_InVar_EIN_0_JuR.mat')
% fixpointbranch=BIF_InVar_EIN_0_JuR{1}
% hopfbranch_psol_1=BIF_InVar_EIN_0_JuR{2}
% hopfbranch_psol_2=BIF_InVar_EIN_0_JuR{3}

hopfbif_branch_pPyEIN_3=df_brnch(funcs,[ind_pPy,ind_pEIN],'hopf'); % use hopf point as first point of hopf branch:

pPy_min= -100;
pPy_max = 350;
pPy_step = 2.5;


pEIN_min = -100;
pEIN_max = 350;
pEIN_step = 1.5;

ind
indexHopf       = 2883%670
hopfpoint=p_tohopf(funcs,fixpointbranch.point(indexHopf));

hopfbif_branch_pPyEIN_3.parameter.min_bound(1:2,:)=[[ind_pPy pPy_min]' [ind_pEIN pEIN_min]']';
hopfbif_branch_pPyEIN_3.parameter.max_bound(1:2,:)=[[ind_pPy pPy_max]' [ind_pEIN pEIN_max]']';
hopfbif_branch_pPyEIN_3.parameter.max_step(1:2,:)=[[ind_pPy pPy_step]' [ind_pEIN pEIN_step]']';
hopfbif_branch_pPyEIN_3.point=hopfpoint;
hopfbif_branch_pPyEIN_3.method.continuation.plot=1;


hopfpoint.parameter(ind_pEIN)=hopfpoint.parameter(ind_pEIN)+0.01e-3; % perturb hopf point
[hopfpoint,success]=p_correc(funcs,hopfpoint,ind_pPy,[],method.point); % correct hopf point, recompute stability
hopfbif_branch_pPyEIN_3.point(2)=hopfpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[hopfbif_branch_pPyEIN_3,s,f,r]=br_contn(funcs,hopfbif_branch_pPyEIN_3,250)           % continue with plotting hopf branch:
hopfbif_branch_pPyEIN_3=br_rvers(hopfbif_branch_pPyEIN_3)                             % reverse Hopf branch
[hopfbif_branch_pPyEIN_3,s,f,r]=br_contn(funcs,hopfbif_branch_pPyEIN_3,250)           % continue in other direction
xlabel('pPy');ylabel('pEIN');

%

bif_hopf_pPypEIN_3 = zeros(2,length(hopfbif_branch_pPyEIN_3.point));
for i=1:length(hopfbif_branch_pPyEIN_3.point)
    bif_hopf_pPypEIN_3(1,i)=hopfbif_branch_pPyEIN_3.point(i).parameter(ind_pPy);
    bif_hopf_pPypEIN_3(2,i)=hopfbif_branch_pPyEIN_3.point(i).parameter(ind_pEIN);
end

figure(2); clf();
plot(bif_hopf_pPypEIN_3(2,:), bif_hopf_pPypEIN_3(1,:))
hold on;

% in the following routine, (a short part of) the periodic solution for a certain point along the
% Hopf-bifurcation branch is calculated. A function fun_determ_Floquet_HpfBifBranch
% is called which provides Floquet-Multipliers. those are used to determine whether the periodic
% solution is stable or not. However, a stable point can rapidly change in
% to a unstable by means of a hopf-fold bifurcation
n=1
for i=61:5:111%length(hopfbif_branch_pPyEIN_3.point)
    i
    floq_hopf3(n,1)=hopfbif_branch_pPyEIN_3.point(i).parameter(ind_pEIN);   % extract current pEIN value
    floq_hopf3(n,2)=hopfbif_branch_pPyEIN_3.point(i).parameter(ind_pPy);   % extract current pPy value
    hopfpoint=hopfbif_branch_pPyEIN_3.point(i);                     % extract current data
    PamRange(1)=hopfpoint.parameter(ind_pEIN) - 0.2*abs(hopfpoint.parameter(ind_pEIN)); % estimate lower border of periodic solution 
    PamRange(3)=hopfpoint.parameter(ind_pEIN) + 0.2*abs(hopfpoint.parameter(ind_pEIN)); % estimate upper border of periodic solution
    PamRange(2)=(abs(PamRange(3))-abs(PamRange(2)))/30;
    temp=fun_determ_Floquet_HpfBifBranch( funcs,hopfpoint,fixpointbranch,indexHopf,ind_pEIN,PamRange,20)
    if length(temp)==2
        floq_hopf3(n,3:4)=temp;
    elseif length(temp)==4
        floq_hopf3(n,3:6)=temp;
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
for i=1:length(floq_hopf3(:,1))
    if max(floq_hopf3(i,3:4))>1.005
        plot(floq_hopf3(i,1),floq_hopf3(i,2),'.r','MarkerSize',20)
    elseif max(floq_hopf3(i,3:4))<1.001    
        plot(floq_hopf3(i,1),floq_hopf3(i,2),'.g','MarkerSize',20)
    else
        plot(floq_hopf3(i,1),floq_hopf3(i,2),'.m','MarkerSize',20)
    end
end





figure(16); 
%clf();
plot(bif_hopf_pPypEIN_3(2,:), bif_hopf_pPypEIN_3(1,:))
plot(bif_hopf_pPypEIN_1(2,:), bif_hopf_pPypEIN_1(1,:))
hold on;
xlabel('pEIN');ylabel('pPy');
legend('hopfbifurcation branch');
title('pPy-pEIN- bifurcation continuation')



% print -dpdf 'InVar_pPypEIN' -r1000 
% bif_InVar_pPypEIN={bif_hopf_pPypEIN_1,floq_hopf1,bif_hopf_pPypEIN_2,bif_hopf_pPypEIN_3,floq_hopf3}
% save('bif_InVar_pPypEIN.mat','bif_InVar_pPypEIN')
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - lower fold - in pEIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_pEIN=18;
ind_pPy=17;


% load('BIF_InVar_EIN_0_JuR.mat')
% fixpointbranch=BIF_InVar_EIN_0_JuR{1}
% hopfbranch_psol_1=BIF_InVar_EIN_0_JuR{2}
% hopfbranch_psol_2=BIF_InVar_EIN_0_JuR{3}

lowerfold_branch_pPypEIN=df_brnch(funcs,[ind_pPy,ind_pEIN],'fold');

pPy_min= -100;
pPy_max = 350;
pPy_step = 2.5;


pEIN_min = -100;
pEIN_max = 350;
pEIN_step = 1.5;

ind
indexLowerFold_pPypEIN       = 824; %howto: check fixpointbranch.point(825).parameter(17) for folds
lowerfoldpoint=p_tofold(funcs,fixpointbranch.point(indexLowerFold_pPypEIN));

lowerfold_branch_pPypEIN.parameter.min_bound(1:2,:)=[[ind_pPy pPy_min]' [ind_pEIN pEIN_min]']';
lowerfold_branch_pPypEIN.parameter.max_bound(1:2,:)=[[ind_pPy pPy_max]' [ind_pEIN pEIN_max]']';
lowerfold_branch_pPypEIN.parameter.max_step(1:2,:)=[[ind_pPy pPy_step]' [ind_pEIN pEIN_step]']';
lowerfold_branch_pPypEIN.point=lowerfoldpoint;
lowerfold_branch_pPypEIN.method.continuation.plot=1;

lowerfoldpoint.parameter(ind_pEIN)=lowerfoldpoint.parameter(ind_pEIN)+0.00001e-3; % perturb hopf point
[lowerfoldpoint,success]=p_correc(funcs,lowerfoldpoint,ind_pPy,[],method.point); % correct hopf point, recompute stability
lowerfold_branch_pPypEIN.point(2)=lowerfoldpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[lowerfold_branch_pPypEIN,s,f,r]=br_contn(funcs,lowerfold_branch_pPypEIN,850)           % continue with plotting hopf branch:
lowerfold_branch_pPypEIN=br_rvers(lowerfold_branch_pPypEIN)                             % reverse Hopf branch
[lowerfold_branch_pPypEIN,s,f,r]=br_contn(funcs,lowerfold_branch_pPypEIN,800)           % continue in other direction
xlabel('pPy');ylabel('pEIN');



bif_lfb_pPypEIN = zeros(2,length(lowerfold_branch_pPypEIN.point));
for i=1:length(lowerfold_branch_pPypEIN.point)
    bif_lfb_pPypEIN(1,i)=lowerfold_branch_pPypEIN.point(i).parameter(ind_pPy);
    bif_lfb_pPypEIN(2,i)=lowerfold_branch_pPypEIN.point(i).parameter(ind_pEIN);
end

figure(3); clf();
plot(bif_lfb_pPypEIN(2,:), bif_lfb_pPypEIN(1,:))


figure(16); clf();
plot(bif_lfb_pPypEIN(2,:), bif_lfb_pPypEIN(1,:))
hold on;
plot(bif_hopf_pPypEIN_3(2,:), bif_hopf_pPypEIN_3(1,:))
plot(bif_hopf_pPypEIN_1(2,:), bif_hopf_pPypEIN_1(1,:))

xlabel('pEIN');ylabel('pPy');
legend('lower fold bifurcation branch','hopfbifurcation branch');
title(['pPy-pEIN- bifurcation continuation, pext=' num2str(pext) 'pps' ])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2 parameter bifurcation - upper fold - in pEIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: identify relevant parameters and ranges through series of codim 1
% bifurcation diagrams

ind_pEIN=18;
ind_pPy=17;


% load('BIF_InVar_EIN_0_JuR.mat')
% fixpointbranch=BIF_InVar_EIN_0_JuR{1}
% hopfbranch_psol_1=BIF_InVar_EIN_0_JuR{2}
% hopfbranch_psol_2=BIF_InVar_EIN_0_JuR{3}


upperfold_branch_pPypEIN=df_brnch(funcs,[ind_pPy,ind_pEIN],'fold'); 

pPy_min= -100;
pPy_max = 350;
pPy_step = 2.5;


pEIN_min = -100;
pEIN_max = 350;
pEIN_step = 1.5;


ind
indexUpperFold       = 1453; %howto: check fixpointbranch.point(825).parameter(17) for folds
upperfoldpoint=p_tofold(funcs,fixpointbranch.point(indexUpperFold));

upperfold_branch_pPypEIN.parameter.min_bound(1:2,:)=[[ind_pPy pPy_min]' [ind_pEIN pEIN_min]']';
upperfold_branch_pPypEIN.parameter.max_bound(1:2,:)=[[ind_pPy pPy_max]' [ind_pEIN pEIN_max]']';
upperfold_branch_pPypEIN.parameter.max_step(1:2,:)=[[ind_pPy pPy_step]' [ind_pEIN pEIN_step]']';
upperfold_branch_pPypEIN.point=upperfoldpoint;
upperfold_branch_pPypEIN.method.continuation.plot=1;

upperfoldpoint.parameter(ind_pEIN)=upperfoldpoint.parameter(ind_pEIN)+0.001e-3; % perturb hopf point
[upperfoldpoint,success]=p_correc(funcs,upperfoldpoint,ind_pPy,[],method.point); % correct hopf point, recompute stability
upperfold_branch_pPypEIN.point(2)=upperfoldpoint;                                 % use as second point of hopf branch:


figure(6); clf;
[upperfold_branch_pPypEIN,s,f,r]=br_contn(funcs,upperfold_branch_pPypEIN,550)           % continue with plotting hopf branch:
upperfold_branch_pPypEIN=br_rvers(upperfold_branch_pPypEIN)                             % reverse Hopf branch
[upperfold_branch_pPypEIN,s,f,r]=br_contn(funcs,upperfold_branch_pPypEIN,1100)           % continue in other direction
xlabel('pPy');ylabel('pEIN');



bif_ufb_pPypEIN = zeros(2,length(upperfold_branch_pPypEIN.point));
for i=1:length(upperfold_branch_pPypEIN.point)
    bif_ufb_pPypEIN(1,i)=upperfold_branch_pPypEIN.point(i).parameter(ind_pPy);
    bif_ufb_pPypEIN(2,i)=upperfold_branch_pPypEIN.point(i).parameter(ind_pEIN);
end

figure(3); clf();
plot(bif_ufb_pPypEIN(2,:), bif_ufb_pPypEIN(1,:))


figure(16); clf();
plot(bif_ufb_pPypEIN(2,:), bif_ufb_pPypEIN(1,:),'-.b')
hold on;
plot(bif_lfb_pPypEIN(2,:), bif_lfb_pPypEIN(1,:))
plot(bif_hopf_pPypEIN_3(2,:), bif_hopf_pPypEIN_3(1,:))
plot(bif_hopf_pPypEIN_1(2,:), bif_hopf_pPypEIN_1(1,:))

xlabel('pEIN');ylabel('pPy');
legend('hopfbifurcation branch','hopfbifurcation branch''upper fold bifurcation branch','lower fold bifurcation branch');
title(['pPy-pEIN- bifurcation continuation' ])






print -dpdf 'InVar_pPypEIN' -r1000 
bif_InVar_pPypEIN={bif_hopf_pPypEIN_1,floq_hopf1,bif_hopf_pPypEIN_2,bif_hopf_pPypEIN_3,floq_hopf3,bif_ufb_pPypEIN,bif_lfb_pPypEIN}
save('bif_InVar_pPypEIN.mat','bif_InVar_pPypEIN')


bif_InVar_pPypEIN=load('bif_InVar_pPypEIN.mat')
bif_hopf_pPypEIN_1=bif_InVar_pPypEIN{1}
floq_hopf1=bif_InVar_pPypEIN{2}
bif_hopf_pPypEIN_2=bif_InVar_pPypEIN{3}
bif_hopf_pPypEIN_3=bif_InVar_pPypEIN{4}
floq_hopf3=bif_InVar_pPypEIN{5}
bif_ufb_pPypEIN=bif_InVar_pPypEIN{6}
bif_lfb_pPypEIN=bif_InVar_pPypEIN{7}


