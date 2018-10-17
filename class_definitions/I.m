% WDF classes

% class for rapresentation of real current intensity generator

classdef I < OnePort
    
    properties
        A; % source current
    end
    
    methods
        function obj = I(A, PortRes) % if A = 0 -------> open-circuit
            obj.A = A;
            obj.PortRes = PortRes;
            obj.WD = 0;
        end
        
        function WU = WaveUp(obj)
            WU = obj.WD - 2*PortRes*obj.A;
            obj.WU = WU;
        end
    end
end