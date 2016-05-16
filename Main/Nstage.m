function p=Nstage(vf,beta,epsilon,alpha);
%Program for solving the multi-stage (N) rocket equation for the first-stage
%payload ratio, p. (Copyright @2006 by Ashish Tewari)
%beta=(Nx1) vector of ratios of specific impulses to that of the first stage
%(beta(k)=Isp_k/Isp_1); first element of beta should be 1.0
%alpha= (Nx1) vector of ratios of payload ratios to that of the first stage
%(alpha(k)=lambda_k/lambda_1); first element of alpha should be 1.0
%epsilon= (Nx1) vector of structural ratios of the stages
%vf=ratio of total velocity impulse to exhaust speed of first stage
%(vf=Delta_v/v_e1)
%(c) 2006 Ashish Tewari
N=size(beta,1);
p=0.1;
f=vf;
tol=1e-9;
for k=1:N
f=f+beta(k)*log(epsilon(k)+alpha(k)*(1-epsilon(k))*p);
end
while abs(f)>tol
f=vf;
fp=0;
for k=1:N
f=f+beta(k)*log(epsilon(k)+alpha(k)*(1-epsilon(k))*p);
fp=fp+alpha(k)*beta(k)/(epsilon(k)+alpha(k)*(1-epsilon(k))*p);
end
d=-f/fp;
p=p+d;
disp(p)
end