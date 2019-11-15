
%------------------------主程序------------------------ 
function sampling_opticfiber(delta_n,T,L_all,Period)
    %-----------------------------------------------  
    % M-----------------取样个数  
    % delta_n-----------调制折射率变化  
    % T-----------------取样光栅占空比  
    % L_all-------------取样光栅总长度  
    % L_grating---------取样光栅节点长度     
    % Period------------取样光栅的取样周期  
    %------------------------------------------------    
    disp('调制折射率变化为：');disp(delta_n);  
    disp('取样光栅占空比：');disp(T);  
    disp('取样光栅总长度:');disp(L_all);  
    disp('取样光栅周期:');disp(Period);    
    lambda=1e-9*linspace(1548,1552,1024);  
    n_eff=1.45;  
    L_grating=Period*T; 
    M=L_all/Period; 
    disp('取样个数：');disp(M);  
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
%-----------------取样光栅正常段--------------------------- 
function [phai]=Normal_OpticFiber(L_normal,n_eff,j,lambda)     
    beita(j)=2*pi/lambda(j)*n_eff;  
    phai=[exp(-i*beita(j)*L_normal),0;0,exp(i*beita(j)*L_normal)]; 
end 
%-----------------------取样光栅均匀段------------------------ 
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