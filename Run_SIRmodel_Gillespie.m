%Script for running and plotting the stochastic SIR model using Gillespie's 
%algorithm 


%Define model parameters as a structure
para = struct('beta',1/3,'gamma',1/7,'N',1000); 

%Define initial conditions as a structure
ICs = struct('S',para.N-1,'I',1,'R',0);

%Define time to run model for
maxtime = 200;

%Run model by calling function ODE_SIRmodel.m
[Classes_ODE] = ODE_SIRmodel(para,ICs,maxtime);

%Run stochastic model using Gillespie
[Classes_G] = Gillespie_SIRmodel(para,ICs,maxtime);

%Store determinisitc simulation metrics

%Peak size
Peak_height_ODE=max(Classes_ODE.I);
%Final size
Final_size_ODE = Classes_ODE.R(end);
%Duration (remember we have to define a proxy for when the epidemic ends)
[,ix]=find(Classes_ODE.I>1,1,'last')
Duration_ODE=Classes_ODE.t(max([ix 1]));

%Store stochastic simulation metrics
%We'll be storing outputs from multiple realisations as a vector

%Peak size
Peak_height(1)=max(Classes_G.I);
%Final size
Final_size(1) = Classes_G.R(end);
%Duration
[,ix]=find(Classes_G.I==0,1,'first')
Duration(1)=Classes_G.t(ix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot dynamics
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',2)

%Make a figure here which plots stochastic infections by time with a red line
%and ODE outputs as dashed black lines
figure(1)
clf
plot(Classes_G.t,Classes_G.I,'r',Classes_ODE.t,Classes_ODE.I,'--k')
xlabel('Time (days)')
ylabel('Infections (I)')
legend('Gillespie','ODE')


%Next we will plot 100 realisations of the stochastic model against the ODE
%model output

figure(2)
clf
%Plot your existing stochastic realisation on the graph now
plot(Classes_G.t,Classes_G.I,'r')
hold on

%Use a for loop to run each simulation, store metrics then plot the
%dynamics
NRuns=1000;
for i=2:NRuns
    [Classes_G] = Gillespie_SIRmodel(para,ICs,maxtime);
    
    %Store stochastic simulation metrics

    %Peak size
    Peak_height(i)=max(Classes_G.I);
    %Final size
    Final_size(i) = Classes_G.R(end);
    %Duration
    [,ix]=find(Classes_G.I==0,1,'first');
    Duration(i)=Classes_G.t(ix);
    
    h1=plot(Classes_G.t,Classes_G.I,'color','r');
end
    
%Add the ODE solution back on top of the stochastic solutions
h2=plot(Classes_ODE.t,Classes_ODE.I,'--k','linewidth',3);
axis([0 200 0 Inf])
xlabel('Time (days)')
ylabel('Infections (I)')

legend([h1 h2],'Gillespie','ODE')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot histograms of: Size of peak, duration and final size

%Size of peak
figure(3)
clf
histogram(Peak_height,'numbins',30,'normalization','probability')
xlabel('Peak height')
ylabel('Probability')
hold on
ymax=get(gca,'ylim');
plot([Peak_height_ODE Peak_height_ODE], [0 ymax(2)])
legend('1000 Gillespie realisations','ODE')

%Size of outbreak
figure(4)
clf
histogram(Final_size,'numbins',50,'normalization','probability')
xlabel('Final size')
ylabel('Probability')
hold on
ymax=get(gca,'ylim');
plot([Final_size_ODE Final_size_ODE], [0 ymax(2)])
legend('1000 Gillespie realisations','ODE')

%Duration
figure(5)
clf
histogram(Duration,'numbins',50,'normalization','probability')
xlabel('Duration')
ylabel('Probability')
hold on
ymax=get(gca,'ylim');
plot([Duration_ODE Duration_ODE], [0 ymax(2)])
legend('1000 Gillespie realisations','ODE')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute analytical probability of no more cases other than index
P_noinf=para.gamma/(para.gamma+para.beta)

%What percentage of the simulations had a final size of 1?
Final_size1=100*sum(Final_size==1)/NRuns


