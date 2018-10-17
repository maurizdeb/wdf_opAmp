% WDF classes

% class for one port capacitor

classdef C < OnePort
    
    properties
        State;
    end
    
    methods 
        function obj = C(PortRes)
            obj.PortRes = PortRes;
            obj.State = 0;
        end
        
        function WU = WaveUp(obj)
            % this is the unit-time delay. Because we update a with the value of b and then we update b 
            % (through the set method in the abstract class) for the time t+1
            WU = obj.State;
            obj.WU = WU;
        end
    end
end