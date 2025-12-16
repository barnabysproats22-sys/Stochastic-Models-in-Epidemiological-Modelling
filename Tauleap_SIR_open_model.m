%Stochastic (tau leap) SIR model with demography code

function [Classes] = Tauleap_SIR_open_model(para,ICs,maxtime,tau)


%Run the tauleap algorithm for an SIR model

%Store the starting point of the simultion from the ICs and copy to a new
%structure called Classes. Define the starting time to be 0
Classes = ICs;
Classes.t = 0;

%Define the current state of the model from the ICs
S=ICs.S;
I=ICs.I;
R=ICs.R;
%Define the current time in the model as 0
t=0;

%Run the model until either the maxtime is exceeded or until there are no
%infected people remaining
while ((t<maxtime) && (I>0))

%Define event rates 
infection = para.beta*I*S/para.N;
recovery = para.gamma*I;
Sdeath = para.mu*S;
Ideath = para.mu*I;
Rdeath = para.mu*R;

%Compute how many events occur for each time step
inf_events=poissrnd(tau*infection);
recovery_events=poissrnd(tau*recovery);
Sdeath_events=poissrnd(tau*Sdeath);
Ideath_events=poissrnd(tau*Ideath);
Rdeath_events=poissrnd(tau*Rdeath);





%Update events

S = S + Sdeath_events + Ideath_events + Rdeath_events - inf_events - Sdeath_events;
I = I + inf_events - recovery_events - Ideath_events;
R = R + recovery_events - Rdeath_events;


%Check nothing less than zero
if S<0 %"undo" infection
        tmp=S;
        S=0;
        I=I+tmp; %take back out of R class
end
if I<0 %"undo" recovery
        tmp=I;
        I=0;
        R=R+tmp; %take back out of R class
end
if R<0 %"undo" death
    tmp=R;
    R=0;
    S=S+tmp; %take back out of S class
end

%Update time
t = t+tau;

%Save information in the Classes structure by extending the matrices of the
%model state and the associated time
Classes.t = [Classes.t t];
Classes.S = [Classes.S S];
Classes.I = [Classes.I I];
Classes.R = [Classes.R R];
end
    
