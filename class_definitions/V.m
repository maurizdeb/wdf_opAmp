% WDF classes

% class for real voltage generator, with PortRes. For ideal PortRes = 0

classdef V < OnePort
    
    properties
        E; % source voltage. Se E = 0 ----> rappresentazione del cortocircuito
    end
    
    methods
        function obj = V(E, PortRes) % se PortRes = 0 -----> generatore ideale
            obj.E = E;
            obj.PortRes = PortRes; 
            obj.WD = 0; % valore iniziale del generatore
        end
        function WU = WaveUp(obj)
            WU = 2*obj.E - obj.WD;
            obj.WU = WU;
        end
    end
    
end