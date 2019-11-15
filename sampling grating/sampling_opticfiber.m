
%------------------------������------------------------ 
function sampling_opticfiber(delta_n,T,L_all,Period)
    %-----------------------------------------------  
    % M-----------------ȡ������  
    % delta_n-----------���������ʱ仯  
    % T-----------------ȡ����դռ�ձ�  
    % L_all-------------ȡ����դ�ܳ���  
    % L_grating---------ȡ����դ�ڵ㳤��     
    % Period------------ȡ����դ��ȡ������  
    %------------------------------------------------    
    disp('���������ʱ仯Ϊ��');disp(delta_n);  
    disp('ȡ����դռ�ձȣ�');disp(T);  
    disp('ȡ����դ�ܳ���:');disp(L_all);  
    disp('ȡ����դ����:');disp(Period);    
    lambda=1e-9*linspace(1548,1552,1024);  
    n_eff=1.45;  
    L_grating=Period*T; 
    M=L_all/Period; 
    disp('ȡ��������');disp(M);  
    lambda_Brag=1550e-9;     
    L_normal=Period-L_grating;  
    for j=1:1024     
        F1=Optic_Fiber(lambda,lambda_Brag,delta_n,n_eff,j,L_grating);   
        phai=Normal_OpticFiber(L_normal,n_eff,j,lambda);      
        F1=(phai*F1*phai)^M; 
        R(j)=(abs(-F1(2,1)/F1(1,1)))^2;    
    end  
    plot(lambda*1e9,R);  
    axis([1548 1552 0 1]);  
    title('The Reflective of Sampling Grating Fiber');  
    xlabel('Wavelength /nm');  
    ylabel('Reflective');     
    axis on;     
    hold on 
end 
%-----------------ȡ����դ������--------------------------- 
function [phai]=Normal_OpticFiber(L_normal,n_eff,j,lambda)     
    beita(j)=2*pi/lambda(j)*n_eff;  
    phai=[exp(-i*beita(j)*L_normal),0;0,exp(i*beita(j)*L_normal)]; 
end 
%-----------------------ȡ����դ���ȶ�------------------------ 
function [F1]=Optic_Fiber(lambda,lambda_Brag,delta_n,n_eff,j,L_grating)  
    delta=1/2 *2*pi*n_eff*(1./lambda-1/lambda_Brag);  
    k=pi./lambda*delta_n;  
    sigma1=2*pi./lambda*delta_n;  
    sigma2=delta+sigma1;  
    s=sqrt(k.^2-delta.^2); 
    f11(j)=(cosh(s(j)*L_grating)-i*sigma2(j)/s(j)*sinh(s(j)*L_grating));     
    f12(j)=-(i*k(j)/s(j)*sinh(s(j)*L_grating));     
    f21(j)=(i*k(j)/s(j)*sinh(s(j)*L_grating)); 
    f22(j)=(cosh(s(j)*L_grating)+i*sigma2(j)/s(j)*sinh(s(j)*L_grating));     
    F1=[f11(j) f12(j);f21(j) f22(j)]; 
end  