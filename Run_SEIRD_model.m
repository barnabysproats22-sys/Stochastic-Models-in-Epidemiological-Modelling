%Run_SEIRD_model
%Script for running programs needed for assignment 2

%IMPORTANT NOTE: Some of the programs take a while to run. Based on my own
%observations I have written scripts which output the approximate time a
%routine may take. My computer is old so this should act as an upper
%bound.

%==============================================||
%Question 1====================================||
%==============================================||

%Q1(b)

%Define model parameters as a structure In units of days
para = struct('beta',5,'sigma',1/10,'gamma',1/3,'pd',0.01,'N',1e5); 

%Define initial conditions as a structure. 90% immune. Take one of 10%
%suseptible to be an infected. 
ICs = struct('S',0.1*para.N - 1,'E',0,'I',1,'R',0.9*para.N,'D',0);

%Define time to run model for
maxtime = 365*2;

%Run model by calling function ODE_SEIRD_model.m
[Classes_ODE] = ODE_SEIRD_model(para,ICs,maxtime);

total_deaths_ODE = round(Classes_ODE.D(end)) - ICs.D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',2)
%Plot dynamics
figure(1)
clf
plot(Classes_ODE.t,Classes_ODE.I)
title('Infections Over Time. Deterministic Model')
xlabel('Time (days)')
ylabel('Number Infected')
xlim([0 maxtime])

%% 
%==============================================||
%Question 2====================================||
%==============================================||

%Q2(b)

%tic
[Classes_G] = Gillespie_SEIRD_model(para,ICs,maxtime); %Same para,maxtime and ICs.
%toc %Can help to estimate run times.

%Plot Dynamics

figure(2)
clf
hold on
plot(Classes_ODE.t,Classes_ODE.I,'--k') %ODE dashed balck
plot(Classes_G.t,Classes_G.I,'-r')      %Stochstic filled red.
hold off
xlim([0 Classes_G.t(end)])
title('Infections Over Time')
xlabel('Time (days)')
ylabel('Number Infected')
legend('Deterministic Model','Stochastic Model Via Gillespie')

%% 

%Q2(c)

%Duration of disease.
duration_days = round(Classes_G.t(end))

%Assume size here means total who were infected.
%These then end up in R and D. Remember to remove ICs.
final_size = Classes_G.R(end) + Classes_G.D(end) - ICs.R - ICs.D

%% 

%Q2(d)

%Set number of runs.
Nruns = 100;
%In previous run found average run time of 0.6 seconds. This gives approx
%runtime in minutes.
approx_run_time_seconds = seconds(Nruns*0.6);
approx_run_time_seconds.Format ='mm:ss';
approx_run_time_seconds

%Initialise vectors to store final size and run time data.
final_size_vector = zeros(Nruns,1); 
run_time = zeros(Nruns,1);

%Fill vectors and plot all runs.
figure(3)
hold on
for i=1:Nruns
    tic;
    [Classes_G] = Gillespie_SEIRD_model(para,ICs,maxtime);
    final_size_vector(i) = Classes_G.R(end) + Classes_G.D(end) - ICs.R - ICs.D;
    run_time(i)=toc;
    plot(Classes_G.t,Classes_G.I,'-r')
end
plot(Classes_ODE.t,Classes_ODE.I,'--k')
hold off
%% 

%Plot histopgrams of data as a subplot
figure(4)
clf

%Histogram of all data.
subplot(3,1,1)
hold on
histogram(final_size_vector,'numbins',100,'normalization','probability')
axis([0 Inf 0 0.7])
%x line where ODE predicts.
xline(Classes_ODE.R(end) + Classes_ODE.D(end) - ICs.R - ICs.D,'r','LineWidth',1.5);
hold off
title('(a) Histogram of Infection Sizes')
xlabel('Size of Epidemic (Individuals)')
ylabel('Probability')
legend('Deterministic Model','Stochastic Model Via Gillespie')

%Histogram of data in range 5000 days plus
subplot(3,1,2)
hold on
histogram(final_size_vector,'numbins',200,'normalization','probability')
axis([5000 Inf 0 0.04])
%x line where ODE predicts.
xline(Classes_ODE.R(end) + Classes_ODE.D(end) - ICs.R - ICs.D,'r','LineWidth',1.5);
hold off
title('(b) Histogram of Infection Sizes 5000+')
xlabel('Size of Epidemic (Individuals)')
ylabel('Probability')
legend('Deterministic Model','Stochastic Model Via Gillespie')

%Histogram of data below 50 days.
subplot(3,1,3)
hold on
histogram(final_size_vector,'numbins',700,'normalization','probability')
axis([0 50 0 0.7])
hold off
title('(c) Histogram of Infection Sizes 0-50')
xlabel('Size of Epidemic (Individuals)')
ylabel('Probability')
legend('Stochastic Model Via Gillespie')

%% 

%Q2(e)

%Set beta to account for reduced contact rate due to inital immunity
beta = 0.1*para.beta;
gamma = para.gamma;
T1 = beta/(beta+gamma);
T2 = gamma/(beta+gamma);

%Predict. See report for details.
PSize_leq3_predicted = T2 + T1*T2^2 + 2*T1^2*T2^3 + 3*T1^3*T2^4
 
%Actual count all with size leq 3.
PSize_leq3_actual=sum(final_size_vector<=3)/Nruns

%In histogram cases that don't explode conatined below 123 so count number
%below this threshold.
P_Extinct_actual = sum(final_size_vector<=123)/Nruns

%Actual is 1/R0, accounting for reduced beta
P_Extinct_predicted = gamma/beta

%% 

%==============================================||
%Question 3====================================||
%==============================================||

%Q3(b)

%Ensuring all aprameters in days^-1.
para2 = struct('beta',5,'sigma',1/10,'gamma',1/3,'mu',1/(60*365),'pd',0.01,'pe',0.95,'v',0,'N',1e6); 
ICs2 = struct('S',66700,'E',430,'I',130,'R',751270,'D',0,'V',181470);maxtime = 20*365; %In days
timestep = 1;
%tic
[Classes_tau] = Tauleap_SEIRDV_open_model(para2,ICs2,maxtime,timestep);
%toc %Again tic toc gives run time allowing prediction of longer runs.
[Classes_SIERDV_ODE] = ODE_SEIRDV_model(para2,ICs2,maxtime);

%Plot one relaisation of tau-leap simulation.
figure(5)
clf
hold on
plot(Classes_tau.t,Classes_tau.I,'r')
plot(Classes_SIERDV_ODE.t,Classes_SIERDV_ODE.I,'--k')
title('Tau Leaping Simulation of Infection Dynamics')
xlabel('Time (Days)')
ylabel('Numebr of People')
legend('Tau Leaping Model','ODE Model')
xlim([0 Classes_tau.t(end)])
hold off
%% 

%Plot multiple realisations using tau-leap
Nruns = 100;

%Approximate total run time. Each run to take on average 7 seconds.
approx_run_time_seconds = seconds(Nruns*7);
approx_run_time_seconds.Format ='mm:ss';
approx_run_time_seconds

%Initialise vectors to store run times and death sizes.
run_time2 = zeros(Nruns,1);
death_sizes = zeros(Nruns,1);

%Plot all runs and populate run_time and death_sizes vectors.
figure(6)
hold on
for i = 1:Nruns
    tic
    [Classes_tau] = Tauleap_SEIRDV_open_model(para2,ICs2,maxtime,timestep);
    death_sizes(i) = Classes_tau.D(end);
    plot(Classes_tau.t,Classes_tau.I,'r')
    run_time2(i) = toc;
end
plot(Classes_SIERDV_ODE.t,Classes_SIERDV_ODE.I,'--k')
hold off
%% 

%Plot death distribution and disease dynamics
figure(7)
clf
set(0,'defaultaxesfontsize',16)
subplot(1,2,1)
plot(Classes_tau.t,Classes_tau.I,'r')
title('Tau-Leaping Simulation of Infection Dynamics');
xlabel('Time (days)')
ylabel('Number Infected')
subplot(1,2,2)
hold on
histogram(death_sizes,'numbins',60,'normalization','probability');
xline(Classes_SIERDV_ODE.D(end),'r','LineWidth',2);
%xlim([2500 4000])
hold off
title('Histogram of Deaths due to Epidemic')
xlabel('Number of Deaths due to Epidemic')
ylabel('Probability')



%% 

%Q3(d)

%Update parameters
para3 = struct('beta',5,'sigma',1/10,'gamma',1/3,'mu',1/(60*365),'pd',0.01,'pe',0.95,'v',0.6,'N',1e6); 

%Set number of runs
Nruns = 300;

%Approximate run time
approx_run_time_seconds = seconds(Nruns*7);
approx_run_time_seconds.Format ='mm:ss';
approx_run_time_seconds

%Initialise vectors to store run times, deaths and duration of runs.
run_time2 = zeros(Nruns,1);
death_sizes = zeros(Nruns,1);
duration = zeros(Nruns,1);

%Plot all simulations. Populate vectors run_time2,death_sizes and duration
figure(8)
clf
hold on
for i = 1:Nruns
    tic
    [Classes_tau] = Tauleap_SEIRDV_open_model(para3,ICs2,maxtime,timestep);
    death_sizes(i) = Classes_tau.D(end);
    duration(i) = Classes_tau.t(end);
    plot(Classes_tau.t,Classes_tau.I,'r')
    run_time2(i) = toc;
end
plot(Classes_SIERDV_ODE.t,Classes_SIERDV_ODE.I,'--k')
hold off

%Count number of runs which end in 5 years.
P5yrs = sum(duration<=5*365)/Nruns

%For calculating v, I simply adapted the parameters in these subroutines.
%My method is described in my report.


%% 

%Figure to plot infection dynamics and histogram of death sizes.
[Classes_tau] = Tauleap_SEIRDV_open_model(para3,ICs2,maxtime,timestep);
[Classes_SEIRDV_ODE] = ODE_SEIRDV_model(para3,ICs2,maxtime);
%% 

figure(9)
clf
set(0,'defaultaxesfontsize',16)
subplot(1,2,1)
hold on
plot(Classes_tau.t,Classes_tau.I,'r')
plot(Classes_SEIRDV_ODE.t,Classes_SEIRDV_ODE.I,'--k')
xlim([0 Classes_tau.t(end)])
title('Tau-Leaping Simulation of Infection Dynamics');
xlabel('Time (days)')
ylabel('Number Infected')
legend('Tau Leaping Simulation', 'ODE Simulation')
hold off
subplot(1,2,2)
hold on
histogram(death_sizes,'numbins',40,'normalization','probability');
xline(Classes_SEIRDV_ODE.D(end),'r','LineWidth',2);
hold off
legend('Tau Leaping Simulation', 'ODE Simulation')
title('Histogram of Deaths due to Epidemic')
xlabel('Number of Deaths due to Epidemic')
ylabel('Probability')


%ReRun Stuff
%Rewrite Answers


