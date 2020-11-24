function x=initialization_m(Nij,xmin,xmax,num_ens)
num_var=3;
num_loc=length(Nij);
x=lhsu(xmin,xmax,num_ens);
x=x';
pop=sum(Nij);
for i=1:num_loc
    x((i-1)*num_var+1,:)=round(x((i-1)*num_var+1,:)*pop(i));%S
    x((i-1)*num_var+2,:)=round(x((i-1)*num_var+2,:)*pop(i));%I
end


function s=lhsu(xmin,xmax,nsample)
% LHS from uniform distribution
nvar=length(xmin);
ran=rand(nsample,nvar);
s=zeros(nsample,nvar);
for j=1: nvar
   idx=randperm(nsample);
   P =(idx'-ran(:,j))/nsample;
   s(:,j) = xmin(j) + P.* (xmax(j)-xmin(j));
end