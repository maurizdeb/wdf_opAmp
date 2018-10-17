% WDF classes

% class for representation of IDEAL current generator

classdef Iideal < WDF
    
    properties 
        A;
        AdaptedRes;
    end
    
    methods
        function obj = Iideal(A, AdaptedRes) %if A = 0 opencircuit
            obj.A = A;
            obj.PortRes = 0;
            obj.AdaptedRes = AdaptedRes;
        end
        
        function WU = WaveUp(obj)
            WU = obj.WD - 2*obj.AdaptedRes*obj.A;
            obj.WU = WU;
        end
    end
    
end