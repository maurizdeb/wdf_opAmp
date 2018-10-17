% WDF classes

% class for linear inductor

classdef L < OnePort
    
    properties
        State;
    end
    
    methods
        function obj = L(PortRes)
            obj.PortRes = PortRes;
            obj.State = 0;
        end
        
        function WU = WaveUp(obj)
            WU = -obj.State;
            obj.WU = WU;
        end
    end
    
end