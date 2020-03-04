
classdef Kalman < handle
    
    properties
        
        nTau;
        nObservation;
        nObservationError = 1e-3;
        
        vKalmanGain;
        vPrediction;
        vEstimate = [0; 0.5];
       
        mA;
        mH = [1, 0];
        mQ = zeros(2);
        mProcessError = 1 * eye(2);
        mStateCovariance;
      
    end
    
    methods
        
        function [obj] = Kalman(nTau)

            obj.nTau = nTau;
            obj.mA = [1, obj.nTau; 0, 1];
            
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