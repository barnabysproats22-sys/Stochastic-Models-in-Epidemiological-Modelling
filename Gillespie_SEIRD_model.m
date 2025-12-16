%Stochastic (Gillespie) SIR model code

function [Classes] = Gillespie_SEIRD_model(para,ICs,maxtime)


%Run Gillespie's algorithm for an SIR model

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
%Define the current time in the model as 0
t=0;

%Run the model until either the maxtime is exceeded or until there are no
%infected people remaining
while ((t<maxtime) && (E+I>0)) %Ends when no I or E that will become I

%Define event rates 
exp = para.beta*I*S/para.N; %Rate move from S to E
inf = para.sigma*E;         %Reate move E to I
rec = (1-para.pd)*para.gamma*I; %Rate inection to recovery
death = para.pd*para.gamma*I;


%Compute the total rate
R_total = exp + inf + rec + death;
%if R_total <=0
%    break
%end

%This conditons seems to prevent everything not starting.
%This over I>0. Unsure how legit but my previous code worked.

%Compute time to next event using a uniform random number
r_1=rand(1);
T_nexttime = -log(r_1)/R_total;

%Select which event has occured
%First infection
r_2=rand(1);

if r_2<exp/R_total
    %Update exposures
    S=S-1;
    E=E+1;
elseif r_2<(exp+rec)/R_total
    %Update recoveries
    I=I-1;
    R=R+1;   
elseif r_2<(exp+rec + death)/R_total
    %Update deaths
    I=I-1;
    D=D+1;
else    %Use elseif so will stop when one is positive
    %Update infections
    E=E-1;
    I=I+1;
end


%Update time
t = t+T_nexttime;

%Save information in the Classes structure by extending the matrices of the
%model state and the associated time
Classes.t = [Classes.t t];
Classes.S = [Classes.S S];
Classes.E = [Classes.E E];
Classes.I = [Classes.I I];
Classes.R = [Classes.R R];
Classes.D = [Classes.D D];
end
    
