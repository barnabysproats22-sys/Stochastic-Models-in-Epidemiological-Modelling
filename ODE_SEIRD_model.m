%ODE SI model code

function [Classes] = ODE_SEIRD_model(para,ICs,maxtime)


%Run ODE using ODE45

%First ODE solver (NB tolerence isn't high enough to see same numerical and
%anlaytical results)
%[t, pop] = ode45(@diff_SIRmodel, [0 maxtime], [ICs.S ICs.I ICs.R], [], para);

%Change ODE options to increase RelTol so that numerical and anlaystical
%results match (to nearest person)
opts = odeset('RelTol',1e-5);
%Only give outputs at defined timepoints
[t, pop] = ode45(@diff_SIRmodel, [0:1: maxtime], [ICs.S ICs.E ICs.I ICs.R ICs.D], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'E', pop(:,2),'I',pop(:,3),'R',pop(:,4),'D',pop(:,5),'t',t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIRmodel(t,pop,para)

%Assign the population matrix into the classes
S=pop(1);
E=pop(2);
I=pop(3);
R=pop(4);
D=pop(5);

%Write down the ODE system
dS = -para.beta*S*I/para.N;
dE = para.beta*S*I/para.N - para.sigma*E;
dI = para.sigma*E - para.gamma*I;
dR = (1-para.pd)*para.gamma*I;
dD = para.pd*para.gamma*I;

%Reshape the derivatives into a column vector
dPop = [dS; dE; dI; dR; dD];

end

end
