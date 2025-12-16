%ODE SI model code

function [Classes] = ODE_SEIRDV_model(para,ICs,maxtime)


%Run ODE using ODE45

%First ODE solver (NB tolerence isn't high enough to see same numerical and
%anlaytical results)
%[t, pop] = ode45(@diff_SIRmodel, [0 maxtime], [ICs.S ICs.I ICs.R], [], para);

%Change ODE options to increase RelTol so that numerical and anlaystical
%results match (to nearest person)
opts = odeset('RelTol',1e-7);
%Only give outputs at defined timepoints
[t, pop] = ode45(@diff_SIRmodel, [0:maxtime], [ICs.S ICs.E ICs.I ICs.R ICs.D ICs.V], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'E', pop(:,2),'I',pop(:,3),'R',pop(:,4),'D',pop(:,5),'V',pop(:,6),'t',t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIRmodel(t,pop,para)

%Assign the population matrix into the classes
S=pop(1);
E=pop(2);
I=pop(3);
R=pop(4);
D=pop(5);
V=pop(6);

%Write down the ODE system
dS = -para.beta*S*I/(para.N-D) + (1-para.pe*para.v)*para.mu*(para.N - D) - para.mu*S;
dE = para.beta*S*I/(para.N-D) - para.sigma*E - para.mu*E;
dI = para.sigma*E - para.gamma*I - para.mu*I;
dR = (1-para.pd)*para.gamma*I - para.mu*R;
dD = para.pd*para.gamma*I;
dV = para.pe*para.v*para.mu*(para.N - D) - para.mu*V;

%Reshape the derivatives into a column vector
dPop = [dS; dE; dI; dR; dD; dV];

end

end
