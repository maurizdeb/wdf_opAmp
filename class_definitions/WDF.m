% WDF classes

% ----------------WDF Class-----------------
classdef WDF < hgsetget
   
    properties
        PortRes % the WDF port resistance
    end
    
    methods
        function Volts = Voltage(obj)
            Volts = (obj.WU + obj.WD)/2; % as defined in the WD literature
        end
    end
    
end