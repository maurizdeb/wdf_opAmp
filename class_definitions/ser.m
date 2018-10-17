% WDF Classes

% Ser Class for bi-series adaptor connection

classdef ser < Adaptor % the class for series 3-port adaptors
   
    properties
        WD % this is the down-going wave at the adapted port
        WU % this is the up-going wave at the adapted port
    end
    
    methods
        function obj = ser(KidLeft, KidRight) % contructor function
            obj.KidLeft = KidLeft; % connect the left child
            obj.KidRight = KidRight; % connect the right child
            obj.PortRes = KidLeft.PortRes + KidRight.PortRes; % adapt. port
        end
        
        function WU = WaveUp(obj) % the up-going wave at the adapted port
            WU = -(WaveUp(obj.KidLeft)+WaveUp(obj.KidRight)); % Wave up
            obj.WU = WU;
        end
        
        function set.WD(obj, WaveFromParent)
            obj.WD = WaveFromParent; % set the down-going wave for the adaptor
            % set the waves to the children according to the scattering
            % rules
            refCoeff = obj.KidLeft.PortRes/obj.PortRes;
            set(obj.KidLeft, 'WD', (obj.KidLeft.WU+obj.KidRight.WU+WaveFromParent)*(-refCoeff)+obj.KidLeft.WU);
            set(obj.KidRight, 'WD', -obj.KidLeft.WU-WaveFromParent + refCoeff*(WaveFromParent + obj.KidLeft.WU+obj.KidRight.WU));
        end
    end
    
end