% WDF classes

% Par Class for bi-parallel adaptor connection

classdef par < Adaptor 
    
    properties
        WU; % for up-going wave
        WD; % for down-going wave
    end
    
    methods
        function obj = par(KidLeft, KidRight)
            obj.KidLeft = KidLeft;
            obj.KidRight = KidRight;
            if obj.KidLeft.PortRes == 0 || obj.KidRight.PortRes == 0
                obj.PortRes = 0;
            else
                obj.PortRes = (1/((1/obj.KidLeft.PortRes) + (1/obj.KidRight.PortRes)));
            end
        end
        
        function WU = WaveUp(obj)
            if obj.KidLeft.PortRes == 0
                refCoeff = 1;
            elseif obj.KidRight.PortRes == 0
                refCoeff = 0;
            else
                refCoeff = (1/obj.KidLeft.PortRes)/(1/obj.PortRes);
            end
            WU = refCoeff*WaveUp(obj.KidLeft) + (1-refCoeff)*WaveUp(obj.KidRight);
            obj.WU = WU;
        end
        
        function set.WD(obj, WaveFromParent)
            obj.WD = WaveFromParent;
            if obj.KidLeft.PortRes == 0
                refCoeff = 1;
            elseif obj.KidRight.PortRes == 0
                refCoeff = 0;
            else
                refCoeff = (1/obj.KidLeft.PortRes)/(1/obj.PortRes);
            end
            set(obj.KidLeft, 'WD', WaveFromParent + (refCoeff-1)*obj.KidLeft.WU+(1-refCoeff)*obj.KidRight.WU);
            set(obj.KidRight, 'WD', WaveFromParent+refCoeff*(obj.KidLeft.WU-obj.KidRight.WU));
        end
    end
    
end