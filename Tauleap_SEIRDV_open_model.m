%Stochastic (tau leap) SIR model with demography code

function [Classes] = Tauleap_SEIRDV_open_model(para,ICs,maxtime,tau)


%Run the tauleap algorithm for an SIR model

%Store the starting point of the simultion from the ICs and copy to a new
%structure called Classes. Define the starting time to be 0
Classes = ICs;
Classes.t = 0;

%Define the current state of the model from the ICs
S=ICs.S;
E=ICs.E;
I=ICs.I;
R=ICs.R;
D=ICs.D;
V=ICs.V;
%Define the current time in the model as 0
t=0;

%Run the model until either the maxtime is exceeded or until there are no
%infected people remaining
while ((t<maxtime) && (E+I>0)) %New condition that no infections can happen

%Define event rates 
new_S = (1-para.pe*para.v)*para.mu*(para.N-D);  %Not constant to ensure total population never changes.
exp = para.beta*I*S/(para.N-D); %Rate move from S to E
inf = para.sigma*E;         %Reate move E to I
rec = (1-para.pd)*para.gamma*I; %Rate inection to recovery
death = para.pd*para.gamma*I;
new_V = para.pe*para.v*para.mu*(para.N-D);

Sdeath = para.mu*S;
Edeath = para.mu*E;
Ideath = para.mu*I;
Rdeath = para.mu*R;
Vdeath = para.mu*V;
%No Ddeath as makes no sense.

%Compute how many events occur for each time step
new_S_events = poissrnd(tau*new_S);
exp_events = poissrnd(tau*exp);
inf_events= poissrnd(tau*inf);
rec_events = poissrnd(tau*rec);
death_events = poissrnd(tau*death);
new_V_events = poissrnd(tau*new_V);

Sdeath_events=poissrnd(tau*Sdeath);
Edeath_events=poissrnd(tau*Edeath);
Ideath_events=poissrnd(tau*Ideath);
Rdeath_events=poissrnd(tau*Rdeath);
Vdeath_events=poissrnd(tau*Vdeath);


%Update events

S = S + new_S_events - Sdeath_events - exp_events;
E = E + exp_events - inf_events - Edeath_events;
I = I + inf_events - Ideath_events - rec_events - death_events;
R = R + rec_events - Rdeath_events;
D = D + death_events;
V = V + new_V_events - Vdeath_events;


%Check nothing less than zero
if S<0 %"undo" exposure
        tmp=S;
        S=0;
        E=E+tmp; %take back out of I class
end
if E<0 %"undo" infection
        tmp=E;
        E=0;
        I=I+tmp; %take back out of I class
end
if I<0 %"undo" recovery
        tmp=I;
        I=0;
        r = rand(1);
        if r<0.99
            R=R+tmp; 
        else
            D=D+tmp;
        end
      %Chance whether we undo a death or recovery. Alt only undo recovery
      %as removing  death is a bit harder to justify.
        
end
if R<0 %"undo" death
    tmp=R;
    R=0;
    S=S+tmp; %take back out of S class
end
if D<0 %"undo" death
    tmp=D;
    D=0;
    S=S+tmp; %take back out of S class
end
if V<0
    tmp = V;
    V=0;
    S = S + tmp; %Instead of losing poeple who are vaccinated, add to those who are not.
end


%Update time
t = t+tau;

%Save information in the Classes structure by extending the matrices of the
%model state and the associated time
Classes.t = [Classes.t t];
Classes.S = [Classes.S S];
Classes.E = [Classes.E E];
Classes.I = [Classes.I I];
Classes.R = [Classes.R R];
Classes.D = [Classes.D D];
Classes.V = [Classes.V V];
end
    
