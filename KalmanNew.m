
classdef KalmanNew < handle
    
    properties
        
        Xk;
        Phi;
        sigma_model;P;
        Q;
        M;
        sigma-Meas;
        R;
        Xk_buffer;
        Z_buffer;
      
    end
    
    methods
        
        function [obj] = KalmanNew(dt)

            Nsamples = length(Xtrue)-1;
    % Previous state (initial guess): Our guess is that the train starts at 0 with velocity
    % that equals to 50% of the real velocity
    obj.Xk = [0; 
    .5*1];
    % Current state estimate
%     Xk=[];
    % Motion equation: Xk = Phi*Xk_prev + Noise, that is Xk(n) = Xk(n-1) + Vk(n-1) * dt
    % Of course, V is not measured, but it is estimated
    % Phi represents the dynamics of the system: it is the motion equation
    obj.Phi = [1 dt;
       0  1];
    % The error matrix (or the confidence matrix): P states whether we should 
    % give more weight to the new measurement or to the model estimate 
    obj.sigma_model = 1;
    % P = sigma^2*G*G';
    obj.P = [obj.sigma_model^2             0;
                 0 obj.sigma_model^2];
    % Q is the process noise covariance. It represents the amount of
    % uncertainty in the model. In our case, we arbitrarily assume that the model is perfect (no
    % acceleration allowed for the train, or in other words - any acceleration is considered to be a noise)
    obj.Q = [0 0;
     0 0];
    % M is the measurement matrix. 
    % We measure X, so M(1) = 1
    % We do not measure V, so M(2)= 0
    obj.M = [1 0];
    % R is the measurement noise covariance. Generally R and sigma_meas can
    % vary between samples. 
    obj.sigma_meas = 1; % 1 m/sec
    obj.R = sigma_meas^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Kalman iteration %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Buffers for later display
    obj.Xk_buffer = zeros(2,Nsamples+1);
    obj.Xk_buffer(:,1) = Xk;
    obj.Z_buffer = zeros(1,Nsamples+1);
            
        end
    
        function [] = iterate(obj, nNewData)
           
            % Observed State
            obj.nObservation = nNewData + obj.nObservationError * randn();

            % State Covariance Update
            obj.mStateCovariance = obj.mA * obj.mProcessError * obj.mA' + obj.mQ;
            
            % Kalman Gain
            obj.vKalmanGain = obj.mStateCovariance * obj.mH' ./ (obj.mH * obj.mStateCovariance * obj.mH' + obj.nObservationError^2);
            
            % Process Error Update
            obj.mProcessError = obj.mStateCovariance - obj.vKalmanGain * obj.mH * obj.mStateCovariance;

            % Estimation
            obj.vEstimate = obj.mA * obj.vEstimate + obj.vKalmanGain * (obj.nObservation - obj.mH * obj.mA * obj.vEstimate); 
                  
        end
        
        function [nTheta, nVelocity] = getData(obj)
            
            nTheta = obj.vEstimate(1);
            nVelocity = obj.vEstimate(2);
            
        end
            
    end
    
end