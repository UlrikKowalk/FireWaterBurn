classdef Tracking < handle
    
    properties
        
        bInitial = true;
        
        nTau;
        mPhi;
        mF;
        
        vG;
        vH;
        
        nManeuverTime = 1;
        nAlpha;
        
        T_constant;
        PP0;
        A;
        rho1;
        rho2;
        rho3;
        B;
        PHI;
        Q_test;
        Q;
        eVL;
        VL;
        C;
        D;

        x0;
        w_std;
        R1;
        Q_sqrt;
        
    end
    
    methods
        
        function [obj] = Tracking(nTau)
            
            obj.nTau = nTau;
            
            obj.nAlpha = 1/obj.nManeuverTime;
            
            obj.mPhi = [1, obj.nTau, 0.5*obj.nTau^2; ...
                0, 1, obj.nTau; ...
                0, 0, 1];
            
            obj.mF = [0, 1, 0; 0, 0, 1; 0, 0, -obj.nAlpha];
            
            obj.vG = [0; 0; 1];
            
            obj.vH = [1, 0, 0];
            
            obj.w_std=0.01*pi/180;      % std of process noise
            obj.R1=(0.1*pi/180)^2;      % variance of measurement
            
            obj.PP0 =  [1 0 0;           % initial estimate of covariance measurement error
                0 ((5*pi/180)/2)^2 0;
                0  0 ((0.5*pi/180))^2];
            
            obj.T_constant=1;
            obj.rho1=1/T_constant;
            obj.rho2=0;
            obj.rho3=0;
            
            obj.A = [-rho3  1     0;  % continuous time model
                0  -rho2  1;
                0  0  -rho1];
            
            obj.B=[0 0 1]'*w_std;
            
            obj.PHI= expm(A*T);
            
            obj.VL=[-A B*B'; zeros(3,3) A']*T;  % Van Loan identity
            
            
            obj.eVL=expm(VL);
            obj.Q=obj.PHI*obj.eVL(1:3,4:6);
            obj.Q=obj.PHI*obj.eVL(1:3,4:6);
            
            obj.Q_test=[obj.nTau^5/20  obj.nTau^4/8 obj.nTau^3/6;  % for large T
                obj.nTau^4/8  obj.nTau^3/3 obj.nTau^2/2;
                obj.nTau^3/6  obj.nTau^2/2    obj.nTau]*obj.w_std^2;
            obj.C=[1  0   0];
            obj.D=0;
            obj.x0 =[0, 0*pi/180 , 0.0*pi/180];
            
            obj.Q_sqrt=chol(obj.Q);

            
        end
        
        function [] = iterate(nNewData)
            
            if obj.bInitial
                PP = obj.PP0;  % predicted in|in-1
                x1_pred=obj.x0';
            else
                PP = obj.PHI*obj.PPU*obj.PHI' + obj.Q;  % P predicted in|in-1
                x1_pred=obj.PHI*x1_hat;     % predict x1
            end
            
            if obj.bInitial
                obj.bInitial = false;
                x1(n,1:3)=x_init;
            else
                x1(n,1:3)=PHI*x1(n-1,1:3)' + Q_sqrt'*randn(3,1);
            end
            y1(n)=C*x1(n,1:3)';
            y(n)=y1(n) + R1^0.5*randn(1,1);
            
            %             if (n>1)
            e1=(PHI*x1(n-1,:)'-x1(n,:)');
            Q_array(:,:,n)=e1*e1';
            %             end
            
            y_pred=C*x1_pred;   % predicted measurement
            innov = y(n)-y_pred; % innovation
            
            S=C*PP*C' + R1;
            K = PP*C'*pinv(S);  % kalman gain
            
            x1_hat = x1_pred + K*innov;  % estimate of x1
            PPU = (eye(3) - K*C)*PP;    % P updated in|in
            
            x1_hat_array(:,n)=x1_hat;
            PPU_array(:,:,n)=PPU;
            gain_array(:,n)=K;
            innov_array(n)=innov;
            
        end
        
        
    end
    
    
    
end