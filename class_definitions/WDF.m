% WDF classes

% ----------------WDF Class-----------------
classdef WDF < hgsetget
   
    properties
        PortRes; % the WDF port resistance
        PortCond;
    end
    
    methods
        function Volts = Voltage(obj)
            Volts = (obj.WU + obj.WD)/2; % as defined in the WD literature
        end
        
        function PortCond = Conductance(obj)
            PortCond = 1/obj.PortRes;
            obj.PortCond = PortCond;
        end
    end
    
end