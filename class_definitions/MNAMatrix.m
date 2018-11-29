classdef MNAMatrix < hgsetget
    %MNAMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X;
        numNodes;
        numPort;
        numAbsorbedElements;
    end
    
    methods
        %constructor
        function obj = MNAMatrix(numNodes, numPort, numAbsorbedElements)
            obj.numNodes = numNodes;
            obj.numPort = numPort;
            obj.numAbsorbedElements = numAbsorbedElements;
            
            dim = numNodes+numPort+numAbsorbedElements;
            
            obj.X = sym(zeros(dim, dim));
        end
        
        function addResistorStamp(obj, Rval, i, j)
            obj.X([i, j], [i, j]) = obj.X([i, j], [i, j]) + [1/Rval, -1/Rval; -1/Rval, 1/Rval];
        end

        function addVoltageSourceStamp(obj, i, j, numSource)
            % As a convention, the + terminal of the voltage source is
            % always i, the - terminal is always j.
            obj.X([i, j, obj.numNodes+numSource], [i, j, obj.numNodes+numSource]) = obj.X([i, j, obj.numNodes+numSource], [i, j, obj.numNodes+numSource]) + [0, 0, 1; 0, 0, -1; 1, -1, 0];
        end  
        
        function addVCVSStamp(obj, i, j, k, l, numSource, gain)
            % As a convention, the + terminal of the control voltage is
            % always i, the - terminal is always j.
            %the + terminal of the controlled voltage source is k, the
            %negative is l
            obj.X([k, l, obj.numNodes+numSource], [i, j, k, l, obj.numNodes+numSource]) = obj.X([k, l, obj.numNodes+numSource], [i, j, k, l, obj.numNodes+numSource]) + [0, 0, 0, 0, 1; 0, 0, 0, 0, -1; -gain, gain, 1, -1, 0];
        end  
        
        function addNullorStamp(obj, i, j, k, l, numSource)
            % As a convention i is the + of the norator, j is the - of the
            % norator
            % k is the + of the nullator, l is the - of the nullator
            obj.X([k, l, i, j, obj.numNodes + numSource],[k, l, i, j, obj.numNodes+numSource]) = obj.X([k, l, i, j, obj.numNodes + numSource],[k, l, i, j, obj.numNodes+numSource]) + [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 1; 0, 0, 0, -1; +1, -1, 0, 0];
        end
        
    end
    
end
