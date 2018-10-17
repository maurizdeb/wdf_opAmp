% WDF classes

% class for OnePort element

classdef OnePort < WDF
    
    properties 
        WU;
        WD;
    end
    
    methods
        function set.WD(obj, val)
            obj.WD = val;
            if or(isa(obj, 'C'), isa(obj, 'L'))
                obj.State = val;
            end
        end
    end
    
end