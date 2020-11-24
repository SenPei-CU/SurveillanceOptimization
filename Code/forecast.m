function forecast(fweek)
%code for Pei et al. Optimizing respiratory virus surveillance networks 
%using uncertainty propagation
%generate 1- to 4-week ahead ILI+ predictions starting from fweek
%fweek: forecast week (fweek <= 27)
%surveillance networks can be modified by changing the variable "observed"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load outbreak data. ilip: ILI+ in 35 states in 31 weeks; ili: ILI in 35
%states in 31 weeks; startweek: the first week; ntest: the number of lab
%tests in 35 states in 31 weeks; pr: test positivity rate in 35 states in
%31 weeks
load Outbreakdata
%load the daily absolute humidity in 35 states
load AH35state
%load the commuting between 35 states
load C_35state
%load scaling parameters for 35 states
load scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set observed states, states are listed in state.mat
observed=1:35;
%set parameters for observational error variance (OEV)
sigma=0.5^2;
nu=0.2^2;
num_loc=size(C,1);%number of locations
tmstep=7;%weekly observation
dt=1;%daily integration
num_times=size(ilip,1);%number of weeks
ts=7*startweek;%starting date
num_ens=300;%number of ensemble members
lambda=1.05;%inflation factor
discrete=0;%run model deterministically
%%%%%%%%%%%%%%%%%%%%initial parameter range
Sl=0.5; Su=0.85;
Il=0.0005; Iu=0.001;
Ll=2*365; Lu=10*365;
Dl=5; Du=7;
R0mxl=2.0; R0mxu=2.5;
R0mnl=1.5; R0mnu=2.0;
thetal=0; thetau=4.5;
xmin=[];xmax=[];
for l=1:num_loc
    xmin=[xmin,Sl,Il,0];
    xmax=[xmax,Su,Iu,0];
end
xmin=[xmin,R0mxl,R0mnl,Ll,Dl,thetal];
xmax=[xmax,R0mxu,R0mnu,Lu,Du,thetau];
%%%%%%%%%%%%%%%%%%%%%%%%%
obs=ilip';%set observation
%set OEV
oev=zeros(num_loc,num_times);
for l=1:num_loc
    for t=1:num_times
        ave=mean(pr(max(1,t-2):t,l));
        ilit=ili(t,l);
        ntestt=max(1,ntest(t,l));
        oev(l,t)=ilit^2/ntestt*(sigma+nu*ave^2);
    end
end
n=length(xmax);%size of state vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% start EAKF of network model
%%%prepare the commute network
Nij=C(:,:,ceil(rand()*1000));
%%%initialization
x=initialization_m(Nij,xmin,xmax,num_ens);
%%%run for one week to initialize system
tcnt=ts-7;
for t=1:1 %1 weeks
    [x,obs_ens]=mixnetworkmodelp(x,tcnt,dt,tmstep,Nij,AH,discrete);
    obs_ens=obs_ens./(scale*ones(1,num_ens));
    tcnt=tcnt+tmstep;
end
%%%% Begin looping through observations
xprior=NaN(n,num_ens,num_times);
xpost=xprior;
obsprior=NaN(num_loc,num_ens,num_times);
obspost=NaN(num_loc,num_ens,num_times);
tcnt=ts;
for tt = 1:fweek%assimilate from week 1
    tt
    %%% inflation of x before assimilation to avoid ensemble collapse
    x=mean(x,2)*ones(1,num_ens)+lambda*(x-mean(x,2)*ones(1,num_ens));
    obs_ens=mean(obs_ens,2)*ones(1,num_ens)+lambda*(obs_ens-mean(obs_ens,2)*ones(1,num_ens));
    %%%%save prior
    xprior(:,:,tt)=x;
    obsprior(:,:,tt)=obs_ens;
    %loop through local observations
    for l=1:num_loc
        if ismember(l,observed)%if observed
        %%%%%  Get the variance of the ensemble
        obs_var = oev(l,tt);
        prior_var = var(obs_ens(l,:));
        post_var = prior_var*obs_var/(prior_var+obs_var);
        if prior_var==0
            post_var=0;
            prior_var=1e-3;
        end
        prior_mean = mean(obs_ens(l,:));
        post_mean = post_var*(prior_mean/prior_var + obs(l,tt)/obs_var);
        %%%% Compute alpha and adjust distribution to conform to posterior moments
        alpha = (obs_var/(obs_var+prior_var)).^0.5;
        dy = post_mean + alpha*(obs_ens(l,:)-prior_mean)-obs_ens(l,:);
        %%%  Getting the covariance of the prior state space and
        %%%  observations  (which could be part of state space, e.g. infections)
        %%%  Loop over each state variable
        rr=zeros(1,size(x,1));
        for j=1:size(x,1)
            A=cov(x(j,:),obs_ens(l,:));
            rr(j)=A(2,1)/prior_var;
        end
        dx=rr'*dy;
        %%%  Get the adjusted ensemble and obs_ens
        x = x + dx;
        obs_ens(l,:)=obs_ens(l,:)+dy;
        obs_ens(l,obs_ens(l,:)<0)=0;
        %%%  Corrections to DA produced aphysicalities
        x = checkbound(x,Nij,xmin,xmax);
        end
    end
    xnew = x;
    %%%%%%%%%%%%%
    xpost(:,:,tt)=xnew;
    obspost(:,:,tt)=obs_ens;
    %%%  Integrate forward one time step
    [x,obs_ens]=mixnetworkmodelp(xnew,tcnt,dt,tmstep,Nij,AH,discrete);
    obs_ens=obs_ens./(scale*ones(1,num_ens));
    tcnt=tcnt+tmstep;
end
%%%%%%%%%%%%%%%%%%%%%
%generate forecast
xpred=NaN(n,num_ens,fweek+4);
obspred=NaN(num_loc,num_ens,fweek+4);
%assign previous observation as true observation
for t=1:fweek
    obspred(:,:,t)=obs(:,t)*ones(1,num_ens);
    xpred(:,:,t)=xpost(:,:,t);
end
%run forward for 4 weeks to make forecast
tpred=ts+(fweek-1)*tmstep;
for t=fweek:fweek+3
    [xpred(:,:,t+1),obspred(:,:,t+1)]=mixnetworkmodelp(xpred(:,:,t),tpred,dt,tmstep,Nij,AH,discrete);
    xpred(:,:,t+1) = checkbound(xpred(:,:,t+1),Nij,xmin,xmax);
    obspred(:,:,t+1)=obspred(:,:,t+1)./(scale*ones(1,num_ens));
    tpred=tpred+tmstep;
end
%plot national average forecast
plot(1:fweek,mean(obs(:,1:fweek),1),'o');
hold on
plot(fweek+1:fweek+4,squeeze(mean(obspred(:,:,fweek+1:fweek+4),1)),'-','Color',[0.8,0.8,0.8]);
plot(fweek+1:fweek+4,mean(obs(:,fweek+1:fweek+4),1),'ro')
xlim([1,fweek+4])
xlabel('Week');
ylabel('ILI+');
%output short-term forecast
%obspred: num_loc x num_ens x fweek+4
save forecastresult.mat obspred
    

function x = checkbound(x,Nij,xmin,xmax)
num_loc=length(Nij);
num_var=3;
N=sum(Nij);
for i=1:num_loc
    %S
    x((i-1)*num_var+1,x((i-1)*num_var+1,:)<0)=1;
    x((i-1)*num_var+1,x((i-1)*num_var+1,:)>N(i))=N(i)-1;
    %I
    x((i-1)*num_var+2,x((i-1)*num_var+2,:)<0)=1;
    x((i-1)*num_var+2,x((i-1)*num_var+2,:)>N(i))=N(i)-1;
    %incidence
    x((i-1)*num_var+3,x((i-1)*num_var+3,:)<0)=1;
    x((i-1)*num_var+3,x((i-1)*num_var+3,:)>N(i))=N(i)-1;
end
%R0max
Rxidx=num_loc*num_var+1;
Rnidx=num_loc*num_var+2;
x(Rxidx,x(Rxidx,:)>xmax(Rxidx))=xmax(Rxidx);
x(Rxidx,x(Rxidx,:)<xmin(Rxidx))=xmin(Rxidx);
%R0min
x(Rnidx,x(Rnidx,:)>xmax(Rnidx))=xmax(Rnidx);
x(Rnidx,x(Rnidx,:)<xmin(Rnidx))=xmin(Rnidx);
%if R0max < R0min
x(Rxidx,x(Rxidx,:)<x(Rnidx,:))=x(Rnidx,x(Rxidx,:)<x(Rnidx,:))+0.01;
%L
Lidx=num_loc*num_var+3;
x(Lidx,x(Lidx,:)>xmax(Lidx))=xmax(Lidx);
x(Lidx,x(Lidx,:)<xmin(Lidx))=xmin(Lidx);
%D
Didx=Lidx+1;
x(Didx,x(Didx,:)>xmax(Didx))=xmax(Didx);
x(Didx,x(Didx,:)<xmin(Didx))=xmin(Didx);
thetaidx=Didx+1;
if num_loc>1
    %q
    x(thetaidx,x(thetaidx,:)<0)=0;
    x(thetaidx,x(thetaidx,:)>xmax(thetaidx))=xmax(thetaidx);
end