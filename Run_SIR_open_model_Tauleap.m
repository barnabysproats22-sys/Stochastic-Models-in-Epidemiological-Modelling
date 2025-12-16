%Script for running and plotting the stochastic SIR model using Gillespie's 
%algorithm 


%Define model parameters as a structure
para = struct('beta',1/365,'gamma',1/(365*5),'mu',1/(365*50),'N',1000); 

%Compute R_0
R_0 = para.beta/(para.gamma + para.mu);

%Endemic equilibrium
S_star = para.N/R_0;
I_star = para.mu*para.N*(R_0 - 1)/para.beta;
R_star = para.N - S_star - I_star;

%Define initial conditions as a structure
ICs = struct('S',S_star,'I',I_star,'R',R_star);

%Update model parameters to account for reducing in infectious duration
para.gamma = 1/365;

%Define time to run model for
maxtime = 20*365;
%Define the time step in tau leaping 
timestep = 1;

%Run model by calling function ODE_SIRmodel.m
[Classes_ODE] = ODE_SIR_open_model(para,ICs,maxtime);

%Run stochastic model using tau leap
tic
[Classes_TL] = Tauleap_SIR_open_model(para,ICs,maxtime,timestep);
toc

%Store metrics
%Time to extinction from ODE (proxy cut of of <1 infection)

    [,ix]=find(Classes_ODE.I>1,1,'last');
    Time_to_ext_yrs_ODE=Classes_ODE.t(max([ix 1]))/365;
    
%Time to extinction from tau leap (stochastic model)

    [,ix]=find(Classes_TL.I==0,1,'first');
    Time_to_ext_yrs(1)=Classes_TL.t(ix)/365;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot dynamics
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',2)

%Make a figure here which plots stochastic infections by time with a red line
%and ODE outputs as dashed black lines
figure(1)
clf
plot(Classes_TL.t./365,Classes_TL.I,'r',Classes_ODE.t./365,Classes_ODE.I,'--k')
xlabel('Time (years since new treatment)')
ylabel('Infections (I)')
legend('Tau leap','ODE')
axis([0 maxtime/365 0 Inf])

%Next we will plot 1000 realisations of the stochastic model against the ODE
%model output
figure(2)
clf
%Plot your existing stochastic realisation on the graph now
plot(Classes_TL.t./365,Classes_TL.I,'r')
hold on

%Use a for loop to run each simulation, store metrics then plot the
%dynamics
NRuns=1000;
tic
for i=2:NRuns
    
    [Classes_TL] = Tauleap_SIR_open_model(para,ICs,maxtime,timestep);
    
    %Store stochastic simulation metrics
    %Time to extinction
    [,ix]=find(Classes_TL.I==0,1,'first');
    Time_to_ext_yrs(i)=Classes_TL.t(ix)/365;
    
    h1=plot(Classes_TL.t./365,Classes_TL.I,'color','r');
   
end
toc
    
%Add the ODE solution back on top of the stochastic solutions
h2=plot(Classes_ODE.t./365,Classes_ODE.I,'--k','linewidth',3);
axis([0 maxtime/365 0 Inf])
xlabel('Time (years since new treatment)')
ylabel('Infections (I)')

legend([h1 h2],'Tau leap','ODE')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot histograms of: Size of peak, duration and final size


%Time to extinction
figure(3)
clf
histogram(Time_to_ext_yrs,'numbins',50,'normalization','probability')
xlabel('Time to extinction (years)')
ylabel('Probability')
hold on
ymax=get(gca,'ylim');
plot([Time_to_ext_yrs_ODE Time_to_ext_yrs_ODE], [0 ymax(2)])
legend('1000 tau leap realisations','ODE')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%What what the median time to extinction and 95% prediction interval?
Median_ext_time = quantile(Time_to_ext_yrs,0.5)

Quant_ext_time = quantile(Time_to_ext_yrs,[0.025 0.975])

ODE_ext_time = Time_to_ext_yrs_ODE
