
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


p_ext_Py=120.0;                   % input to Py 0->no input to Py
p_ext_EIN=-20.0;                   % input to EIN 0-> no input to EIN
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
 

% Delays
neuron_tau=@()[25];

% Bifurcation parameter

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
fixpointbranch=df_brnch(funcs,ind_pEIN,'stst')

% set bounds for continuation parameter
freepammin= -120;
freepammax = 350;
freepamstep = 0.25; %1.5

fixpointbranch.parameter.min_bound(1,:)=[ind_pEIN freepammin];
fixpointbranch.parameter.max_bound(1,:)=[ind_pEIN freepammax];
fixpointbranch.parameter.max_step(1,:)=[ind_pEIN freepamstep];
% use stst as a first branch point:
fixpointbranch.point=stst;


%  Extend and continue branch of trivial equilibria
stst.parameter(ind_pEIN)=stst.parameter(ind_pEIN)+0.001e-3;
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
title(['point: ' num2str(ind_point) ' | parameter:' num2str(fixpointbranch.point(ind_point).parameter(ind_pEIN))])


% plot stability versus point number:
figure(1); clf;
br_plot(fixpointbranch,[],ym,'b');
br_plot(fixpointbranch,[],ym,'b.');
plot([0 30],[0 0],'-.');
xlabel('point number along branch');ylabel('\Re(\lambda)');

% Plot all Eigenvalues along branch and along free parameter

branch_summar=zeros(22,length(fixpointbranch.point));
for i=1:length(fixpointbranch.point)
    branch_summar(1,i)=fixpointbranch.point(i).parameter(ind_pEIN);   % parameter value
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

freePamRange(1) = -105;
freePamRange(2) = 1.5;
freePamRange(3) = 150;
indexHopf       = 1176;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_1 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pEIN,freePamRange,150,1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Continuation of second Hopf bifurcation in a periodic solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind

freePamRange(1) = -108;
freePamRange(2) = 2.5;
freePamRange(3) = 350;
indexHopf       = 1587;       %get from eigenvalues plot or print ind and points

[ hopfbranch_psol_2 ] = fun_cont_Hopfbr( funcs,fixpointbranch,indexHopf,ind_pEIN,freePamRange,150,1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Put Branches together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first hopf branch
[ extval_Hopf_1 ] = fun_prepareBranchPlot( hopfbranch_psol_1 , ind_pEIN);
[ extval_Hopf_2 ] = fun_prepareBranchPlot( hopfbranch_psol_2 , ind_pEIN);



figure(12);clf();
%plot(fix_pot(1,:),fix_pot(2,:)./560,'b')
a=fun_plotFPcurve(fixpointbranch,ind,ind_pEIN,12);
hold on;


plot(extval_Hopf_1(1,:),extval_Hopf_1(2,:)./560,'r');
plot(extval_Hopf_1(1,:),extval_Hopf_1(3,:)./560,'r');

plot(extval_Hopf_2(1,:),extval_Hopf_2(2,:)./560,'g');
plot(extval_Hopf_2(1,:),extval_Hopf_2(3,:)./560,'g');

xlabel('pext EIN')
ylabel('pyrapot')
title(['pext Py:' num2str(p_ext_Py) 's-1 '])



print -dpdf 'BIF_InVar_EINsweep_Py120' -r1000 

BIF_InVar_EINsweep_Py120={fixpointbranch}%,hopfbranch_psol_1,hopfbranch_psol_2}
save('BIF_InVar_EINsweep_Py120.mat','BIF_InVar_EINsweep_Py120')






