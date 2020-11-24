function [x,obs]=mixnetworkmodelp(x,ts,dt,tmstep,Nij,AH,discrete)
%x:state vector,num_var:number of variables for each location,num_loc:
%number of locations,ts:integration start time,dt:minimal integration
%timestep,tmstep:total integration timestep
%Nij(i,j): i<-j, live in j and work at i
%discrete=0: continuous, 1: stochastic
Nave=(Nij+Nij')/2;%average commuting flux
Nave(:,6)=Nij(:,6);
Nave(6,:)=Nij(:,6)';
Cave=zeros(length(Nij),1);%average commuting in each location
for i=1:length(Nij)
    Cave(i)=sum(Nave(:,i))-Nave(i,i);
end
%continuous
num_var=3;
num_loc=length(Nij);
num_ens=size(x,2);
N=sum(Nij);%population
AH(length(AH)+1:2*length(AH),:)=AH;
%R0max
Rxidx=num_loc*num_var+1;
%R0min
Rnidx=num_loc*num_var+2;
%L
Lidx=num_loc*num_var+3;
%D
Didx=num_loc*num_var+4;
%q
thetaidx=num_loc*num_var+5;
%prepare BETA
BT1=zeros(length(AH),num_ens,num_loc);
BETA=zeros(length(AH),num_ens,num_loc);
for i=1:num_loc
    for j=1:num_ens
        Rx=x(Rxidx,j);Rn=x(Rnidx,j);D=x(Didx,j);
        b=log(Rx-Rn); a=-180;
        BT1(:,j,i)=exp(a*AH(:,i)+b)+Rn;
        BETA(:,j,i)=BT1(:,j,i)/D;
    end
end
%transform format from state space vector to matrix
I=zeros(num_loc,num_ens,abs(tmstep)+1);
S=zeros(num_loc,num_ens,abs(tmstep)+1);
Incidence=zeros(num_loc,num_ens,abs(tmstep)+1);
obs=zeros(num_loc,num_ens);
for i=1:num_loc
    I(i,:,1)=x((i-1)*num_var+2,:);
    S(i,:,1)=x((i-1)*num_var+1,:);
end
L=x(Lidx,:);
D=x(Didx,:);
theta=x(thetaidx,:);
sk1=zeros(num_loc,num_ens);
ik1=sk1;ik1i=sk1;
sk2=sk1;ik2=sk1;ik2i=sk1;
sk3=sk1;ik3=sk1;ik3i=sk1;
sk4=sk1;ik4=sk1;ik4i=sk1;
%start integration
tcnt=0;
for t=ts+dt:dt:ts+tmstep
    tcnt=tcnt+1;
    %first step
    for i=1:num_loc
        Eimmloss=dt*((N(i)*ones(1,num_ens)-S(i,:,tcnt)-I(i,:,tcnt))./L);
        Einf=min(dt*(BETA(t,:,i).*S(i,:,tcnt).*I(i,:,tcnt)/N(i)),S(i,:,tcnt));
        Erecov=min(dt*(I(i,:,tcnt)./D),I(i,:,tcnt));
        Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
        
        EIleft=min(dt*theta.*I(i,:,tcnt)*Cave(i)/N(i),I(i,:,tcnt));
        ESleft=min(dt*theta.*S(i,:,tcnt)*Cave(i)/N(i),S(i,:,tcnt));

        EIenter=zeros(1,num_ens);
        ESenter=zeros(1,num_ens);
        for j=1:num_loc
            if j~=i
                EIenter=EIenter+ min(dt*theta.*I(j,:,tcnt)*Nave(i,j)/N(j),I(j,:,tcnt));
                ESenter=ESenter+ min(dt*theta.*S(j,:,tcnt)*Nave(i,j)/N(j),S(j,:,tcnt));
            end
        end
        if discrete==0
            sk1(i,:)=Eimmloss-Einf+ESenter-ESleft;
            ik1(i,:)=Einf-Erecov+EIenter-EIleft;
            ik1i(i,:)=Einf;
        end
        if discrete==1
            l=poissrnd([Eimmloss;Einf;Erecov;ESenter;ESleft;EIenter;EIleft]);
            sk1(i,:)=l(1,:)-l(2,:)+l(4,:)-l(5,:);
            ik1(i,:)=l(2,:)-l(3,:)+l(6,:)-l(7,:);
            ik1i(i,:)=l(2,:);
        end
        
    end
    %second step
    Ts1=S(:,:,tcnt)+round(sk1/2);
    Ti1=I(:,:,tcnt)+round(ik1/2);
    
    for i=1:num_loc
        Eimmloss=dt*((N(i)*ones(1,num_ens)-Ts1(i,:)-Ti1(i,:))./L);
        Einf=min(dt*(BETA(t,:,i).*Ts1(i,:).*Ti1(i,:)/N(i)),Ts1(i,:));
        Erecov=min(dt*(Ti1(i,:)./D),Ti1(i,:));
        Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
        
        EIleft=min(dt*theta.*Ti1(i,:)*Cave(i)/N(i),Ti1(i,:));
        ESleft=min(dt*theta.*Ts1(i,:)*Cave(i)/N(i),Ts1(i,:));

        EIenter=zeros(1,num_ens);
        ESenter=zeros(1,num_ens);
        for j=1:num_loc
            if j~=i
                EIenter=EIenter+ min(dt*theta.*Ti1(j,:)*Nave(i,j)/N(j),Ti1(j,:));
                ESenter=ESenter+ min(dt*theta.*Ts1(j,:)*Nave(i,j)/N(j),Ts1(j,:));
            end
        end
        if discrete==0
            sk2(i,:)=Eimmloss-Einf+ESenter-ESleft;
            ik2(i,:)=Einf-Erecov+EIenter-EIleft;
            ik2i(i,:)=Einf;
        end
        if discrete==1
            l=poissrnd([Eimmloss;Einf;Erecov;ESenter;ESleft;EIenter;EIleft]);
            sk2(i,:)=l(1,:)-l(2,:)+l(4,:)-l(5,:);
            ik2(i,:)=l(2,:)-l(3,:)+l(6,:)-l(7,:);
            ik2i(i,:)=l(2,:);
        end
        
    end
    %third step
    Ts2=S(:,:,tcnt)+round(sk2/2);
    Ti2=I(:,:,tcnt)+round(ik2/2);
    
    for i=1:num_loc
        Eimmloss=dt*((N(i)*ones(1,num_ens)-Ts2(i,:)-Ti2(i,:))./L);
        Einf=min(dt*(BETA(t,:,i).*Ts2(i,:).*Ti2(i,:)/N(i)),Ts2(i,:));
        Erecov=min(dt*(Ti2(i,:)./D),Ti2(i,:));
        Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
        
        EIleft=min(dt*theta.*Ti2(i,:)*Cave(i)/N(i),Ti2(i,:));
        ESleft=min(dt*theta.*Ts2(i,:)*Cave(i)/N(i),Ts2(i,:));

        EIenter=zeros(1,num_ens);
        ESenter=zeros(1,num_ens);
        for j=1:num_loc
            if j~=i
                EIenter=EIenter+ min(dt*theta.*Ti2(j,:)*Nave(i,j)/N(j),Ti2(j,:));
                ESenter=ESenter+ min(dt*theta.*Ts2(j,:)*Nave(i,j)/N(j),Ts2(j,:));
            end
        end
        if discrete==0
            sk3(i,:)=Eimmloss-Einf+ESenter-ESleft;
            ik3(i,:)=Einf-Erecov+EIenter-EIleft;
            ik3i(i,:)=Einf;
        end
        if discrete==1
            l=poissrnd([Eimmloss;Einf;Erecov;ESenter;ESleft;EIenter;EIleft]);
            sk3(i,:)=l(1,:)-l(2,:)+l(4,:)-l(5,:);
            ik3(i,:)=l(2,:)-l(3,:)+l(6,:)-l(7,:);
            ik3i(i,:)=l(2,:);
        end
        
    end

    %fourth step
    Ts3=S(:,:,tcnt)+round(sk3);
    Ti3=I(:,:,tcnt)+round(ik3);
    
    for i=1:num_loc
        Eimmloss=dt*((N(i)*ones(1,num_ens)-Ts3(i,:)-Ti3(i,:))./L);
        Einf=min(dt*(BETA(t,:,i).*Ts3(i,:).*Ti3(i,:)/N(i)),Ts3(i,:));
        Erecov=min(dt*(Ti3(i,:)./D),Ti3(i,:));
        Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
        
        EIleft=min(dt*theta.*Ti3(i,:)*Cave(i)/N(i),Ti3(i,:));
        ESleft=min(dt*theta.*Ts3(i,:)*Cave(i)/N(i),Ts3(i,:));

        EIenter=zeros(1,num_ens);
        ESenter=zeros(1,num_ens);
        for j=1:num_loc
            if j~=i
                EIenter=EIenter+ min(dt*theta.*Ti3(j,:)*Nave(i,j)/N(j),Ti3(j,:));
                ESenter=ESenter+ min(dt*theta.*Ts3(j,:)*Nave(i,j)/N(j),Ts3(j,:));
            end
        end
        if discrete==0
            sk4(i,:)=Eimmloss-Einf+ESenter-ESleft;
            ik4(i,:)=Einf-Erecov+EIenter-EIleft;
            ik4i(i,:)=Einf;
        end
        if discrete==1
            l=poissrnd([Eimmloss;Einf;Erecov;ESenter;ESleft;EIenter;EIleft]);
            sk4(i,:)=l(1,:)-l(2,:)+l(4,:)-l(5,:);
            ik4(i,:)=l(2,:)-l(3,:)+l(6,:)-l(7,:);
            ik4i(i,:)=l(2,:);
        end
        
    end

    S(:,:,tcnt+1)=S(:,:,tcnt)+round(sk1/6+sk2/3+sk3/3+sk4/6);
    I(:,:,tcnt+1)=I(:,:,tcnt)+round(ik1/6+ik2/3+ik3/3+ik4/6);
    Incidence(:,:,tcnt+1)=round(ik1i/6+ik2i/3+ik3i/3+ik4i/6);
    
    
    %observation
    for i=1:num_loc
        obs(i,:)=obs(i,:)+Incidence(i,:,tcnt+1)/N(i);
    end
end

for i=1:num_loc
    %S
    x((i-1)*num_var+1,:)=S(i,:,tcnt+1);
    %I
    x((i-1)*num_var+2,:)=I(i,:,tcnt+1);
    %Incidence
    x((i-1)*num_var+3,:)=obs(i,:);
end