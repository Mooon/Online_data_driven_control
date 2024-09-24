clear all  
%close all 
clc 
global A1 A2 B2 B1 m n T  
%% aircraft engine 
Ts = 0.1; 
 
A1c =[-0.2296, 0.9931; 0.02436, -0.2046]; 
B1c = [-0.0434 -0.01145; -1.73 -0.517];  
[A1, B1] = c2d(A1c,B1c, Ts) 
eigs(A1) 
 
A2c = [-1.175 0.9871; -8.458 -0.8776]; 
B2c = [-0.194, -0.03593; -19.29, -3.803];  
[A2, B2] = c2d(A2c,B2c, Ts) 
eigs(A2) 
 
A_sw{1}=A1; A_sw{2}=A2;  
B_sw{1}=B1; B_sw{2}=B2; 
 
%% data setup and collection 
rng('default') 
m = size(B1,2); 
n= size(A1,1); 
N = (m+1)*n + m; 
T = 2*N -1 
 
magx=0.5; 
x = -magx + (magx+magx).*rand(n,1);                
U = -magx + (magx+magx).*rand(m,T);                 
X = []; 
X = x; 
for i =1:T 
    x = A1*x+B1*U(:,i); 
    X = [X, x]; 
end 


%% after i got the dataset offline U_-1, X_-1, X_0, now go online from k=0 
time=150; 
time_faults = [10; 25;50;67;95;110]; %95 
sw_seq = [2,1,2,1,2,1]; 

As = A1; Bs = B1; eps= 0.001; N=[]; 
sw_con=1; count=1; 
Xcl=[]; Ucl=[]; y_faults=[]; 
x0=x; 
K = U(:,end); 
for t=0:time 
   disp(t) 
   if t==time_faults(sw_con) 
        As=A_sw{sw_seq(sw_con)}; 
        Bs=B_sw{sw_seq(sw_con)}; 
        if sw_con== size(time_faults,1) 
             
        else 
            sw_con=sw_con+1; 
        end 
   end 
   X0 = X(:,end-T:end-1); 
   U0 = U(:,end-T+1:end); 
   X1 = X(:,end-T+1:end);  
   old_K = K; 
   %at time t = now calcolo la K(t) 
   cvx_begin sdp quiet 
    variable gam(1,1) 
    variable P(n,n) symmetric 
    variable Q(T,n) 
    variable L(m,m) symmetric 
    minimize gam 
    subject to 
        [L, U0*Q; Q'*U0', P] >= 0 
        [eye(n)-P, (X1)*Q; Q'*(X1)', -P] <= 0 
         X0*Q == P 
         trace(P)+ trace(L) <= gam 
    cvx_end 
    K = U0*Q/P; 
    if isnan(K) 
        K = old_K; 
    end 
    %construct my u(t) and save u(t), x(t) 
    e = -eps + (eps+eps).*(randn(1,m))'; 
    %E = [E, e]; 
    u = K*x + norm(x)*e; 
        %salvo da k=0 
        Xcl = [Xcl, x];  %<-- x(t) 
        Ucl = [Ucl, u]; 
        N = [N, norm(x)]; 
        %%%% 
    x = As*x+Bs*u;  
    %update x(t+1) and u(t) for the next step 
    U = [U, u]; 
    X = [X, x]; 
end 
 
%% 
colors = {'#ff9100','#1e96fc' }; 
figure 
subplot 211 
box on 
for j = 1:n 
    plot(0:Ts:(size(Xcl,2)-1)*Ts, Xcl(j,:), '.-', 'MarkerSize', 7,'linewidth',1.2, 'color', colors{j}); 
    ylabel('$x$', 'interpreter','latex', 'fontsize',16) 
    %xlabel('time[s]', 'interpreter','latex', 'fontsize',13) 
    hold on 
    grid on 
end 
for i=1:size(time_faults,1) 
    index = time_faults(i); 
    for j=1:n 
        plot((index)*Ts,Xcl(:,index+1), '*','MarkerSize', 11, 'linewidth',1.1, 'color', "#db2b39" ); 
        hold on 
    end 
    legend('$x_1$','$x_2$', 'switches', 'interpreter','latex', 'fontsize',13) 
end 
%% 
figure 
h = subplot(2,1,1); 
coloru={'#487bea', '#da627d'}; 
for j = 1:m 
    plot(0:Ts:(size(Ucl,2)-1)*Ts, Ucl(j,:),'.-', 'MarkerSize', 7,'linewidth',1.2, 'color', coloru{j}); 
    ylabel('$u$', 'interpreter','latex', 'fontsize',16) 
    xlabel('time[s]', 'interpreter','latex', 'fontsize',14) 
    legend('$u_1$','$u_2$', 'interpreter','latex', 'fontsize',13) 
    hold on 
    grid on 
end 
 
%save('aircraft_case')
