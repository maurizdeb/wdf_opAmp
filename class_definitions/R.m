% WDF classes

% class for the resistor

classdef R < OnePort
    
    methods
        function obj = R(PortRes)
            obj.PortRes = PortRes;
        end
        function WU = WaveUp(obj)
            WU = 0;
            obj.WU = WU;
        end
    end
    
end