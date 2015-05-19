function genscenarios
clear all; clear mex; clc;
global rd dd rc dc rh dh 
data_params_single_house
    [rd,dd]=scenarios(pon_dish,pof_dish,Nsamples);
    [rc,dc]=scenarios(pon_clothes,pof_clothes,Nsamples);
    [rh,dh]=scenarios(pon_heater,pof_heater,Nsamples);
% ecost=calcost(2,rd,dd,Pd);
    tau0=[0;0;0];
    A=[];
    b=[];
    Aeq=[];
    beq=[];
    lb=[0;0;0];
    ub=[10;10;10];
    options=[];
    vt=[1;1;1];
  [x,h,exitflag,output,lambda,grad,hessian] = fminconint(@calcost_home,tau0,A,b,Aeq,beq,lb,ub,@nonlcon,options,vt,100);
end


  

% %r1=min(r);
% %r2=max(r);
% %range=r1:r2;
% %Nr=hist(r,range);
% %bar(range,Nr./numel(r));
% %d1=min(r);
% %d2=max(r);
% %range=d1:d2;
% %Nd=hist(duration,range);
% %bar(range,Nd./numel(r));
% % tau=0:10;
% % for i=1:length(tau)
% %     ecost_home(i)=calcost_home(tau(i));
% % end
%
% %plot(tau,ecost_home);
% keyboard
% end
%
%
%
%
 function [e_total,g]=calcost_home(tau)
 data_params_single_house;
 global rd dd rc dc rh dh 
       for i =1:Napp
            ecost_dish = calcost(tau(1),rd,dd,Pd);
            ecost_clothes=calcost(tau(2),rc,dc,Pc);
            ecost_heater=calcost(tau(3),rh,dh,Ph);
       end
       e_total=ecost_dish;+ecost_clothes+ecost_heater;
       tau;
       e_total;
       g=0;
 end
  
function [c,ceq] = nonlcon(tau)
c = [];
ceq = [];
end

function ecost=calcost(tau,r,d,P)
data_params_single_house;
%P=P*ones(size(e_t,1),1);
%P=[P;P;P];
%e_t=[e_t;e_t;e_t];
for i=1:length(r)
    %z1=r(i)+tau+1;
    %z2=z1+d(i)-1;
    t1(i)=r(i)+tau-1;
    t2(i)=t1(i)+d(i);
x=t1(i):0.01:t2(i);
for j=1:length(x)
    y(j)=electricost(x(j));
end  
c_elec(i)=trapz(x,y)*(P/d(i));   
y=0;
if tau > m_n
    c_delay=(m_n*c_n1);
    c_add_delay=c_delay+c_n2*(tau-m_n);
else
    c_add_delay=0;
    c_delay=tau*c_n1;
end
C_T(i)=c_elec(i)+c_delay+c_add_delay;
end
ecost=sum(C_T)/Nsamples;
end

function c_elec=electricost(t)
if t>24 && t<=48
    t=t-24;
elseif t>48 && t<=72
    t=t-48;
end
if t < 10 && t >=0
    c_elec=0.014;
elseif t >=10 && t <22
    c_elec =0.21;
elseif t>=22 && t<=24
    c_elec=0.014;
end
end


function [r,duration]=scenarios(p_on,p_off,Nsamples)
%Code written by Arun Srikanth for OPR project.
%This function returns the request time slot (not time) and duration for an applicance
%based on on-off probabilites.
%Further more two assumptions are made while generating the scenarios
%1. The duration of an equipment does not exceed 24 hours
%2. An equipment is allowed to be used only once.
%The input arguments for the functions are
%p_on  %off to on probability
%p_off %on to off porbability
%No of samples you want to generate
k=1; % Sample no: set it initially to 1.
ton_acc=zeros(Nsamples,1); %Initialize
toff_acc=zeros(Nsamples,1);%Initialize
duration=zeros(Nsamples,1); %Initialize
reject=0;
dcounter=0;% This is the day counter in case one wants to see how many days an equipment is requested,
%however one has to comment out certain parts of the code in order to see this
%due to assumption 1.
t=1;%Set the current clock to 1st time slot, t is the time ticker
while k < Nsamples+1 %Check whether the sample number is less thank k
       if t<=24 %Check whether the timeslot is less than 24
        reject=0;   
        p=rand(1); % choose a random number between 0 and 1
        %---------------------Check whether the equipment is on based on given probability-------------------------------------------------------------------
        if  p_on(t) > p %Check whether the machine is on the time slot t
            ton_acc(k)=t; % If yes accept t
            j=1; %Time ticker which keeps  track of the duration
            count=0;
            toff=ton_acc(k)+1; %increase t by t+1
            %---------------This part of the code finds the duration once the request is made--------------------------------------------------------
            %---------------------------------------------------------------------------------------
            while count < 1 %check whether the count is less than one
                if toff > 24 % if the next time slot happens to be on the next day
                    toff=toff-24; %  cycle
                    dcounter=dcounter+1; %Increase the day counter
                end
                p=rand(1); % Choose a random number between 0 and 1
                if p_off(toff) > p %Check whether the machine is off on time slot t+1
                    toff_acc(k)=toff;%accecpt the time slot if yes
                    duration(k)=j;%collect the duration as a vector
                    t=24; %Reset the clock to the end of the day
                    count=count+1; %Increase the counter to exit the while loop since the equipment is not switched twice on same day
                    day(k)=dcounter;% Keep the day counter
                    dcounter=0;%Set the day counter equal to zero 
                else
                    %----------------Remove samples if the duration is more than a day-----------------------------------------------------------------
                    if toff==ton_acc(k)
                        %j %Uncomment this line to see whether the duration is
                        %24 hours
                        count=count+1; %Exit the while loop
                        t=24;%Reset the clock to the end of the day
                        dcounter=0; %Set the day counter back to zero
                        reject=1;
                    end
                    %Otherwise go to next time slot to see if its off
                    j=j+1;
                    toff=toff+1;
                end
            end
            %--------------------------------------------------------------------------------------------------------------------
            %---------------------End of the part that finds the duration---------------------------------------------
        end
        %-----------------------------------------------------------------------------------------------------------------------
        t=t+1;%Increase the ticker if no equipment is on at time t
    elseif t>24 %if the t exceeds the day
        if reject==1 
            k=k; %If the sample is rejected do not upadate sample number
            t=1; %Restart the clock
        elseif reject==0 % If there is no rejection
             k=k+1; % Collect samples and
            t=1; %reset the clock
        end
    end
end
r=ton_acc';
day=day';
ton_acc=ton_acc';
toff_acc=toff_acc';
%Importance sampling 
Scenario= [r' duration];
Scenario(any(Scenario==0,2),:)=[]; % Elminate samples that dont contribute to the average
r= Scenario(:,1);
duration=Scenario(:,2);
end








