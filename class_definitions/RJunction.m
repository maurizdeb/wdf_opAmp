classdef RJunction < WDF
    %RJUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        S;
        WU;
        WD;
        ConnectedPorts;
        AdaptedIndex;
        WUPs_utn;
    end
    
    methods
        function obj = RJunction(MnaMatrix, R_vect, Ports, adaptedPort)
            %MnaMatrix is a MNAMatrix object
            %
            %R_vect is a symbolic vector, with the port resistances, this
            %must be coherent with the matrix X.
            %
            %R_value is the vector with the known (non adapted) port
            %resistances.
            
            X = MnaMatrix.X;
            numNodes = MnaMatrix.numNodes;
            numPort = MnaMatrix.numPort;
            numAbsorbedElements = MnaMatrix.numAbsorbedElements;
            
            X= X(2:end, 2:end);
            X_inv = inv(X);
            
            %compute scattering matrix: S = I + 2[0 R 0] * X^(-1) * [0 I 0]'
            %COULD NOT WORK WITH NULLORS
            R_diag = diag(R_vect);
            left_term = 2*[zeros(numPort, (numNodes-1) ), R_diag, zeros(numPort, numAbsorbedElements)];
            right_term = [zeros(numPort, (numNodes-1) ), eye(numPort), zeros(numPort, numAbsorbedElements)];
            S = eye(numPort) + left_term*(X_inv*(right_term'));
            
            el = find(R_vect==adaptedPort);
            obj.AdaptedIndex = el;
            R_value = [Ports.PortRes];
            [obj.PortRes, param, cond] = solve(S(el,el)==0, adaptedPort, 'ReturnConditions', true);
            assume(cond);
            obj.PortRes = double(subs(obj.PortRes, R_vect([1:(el-1), (el+1):end]), R_value));
            R_value = [R_value(1:el-1), obj.PortRes, R_value(el:end)];
            obj.S = double(subs(S, R_vect, R_value));
            obj.S(el,el) = 0;
            obj.ConnectedPorts = Ports;
        end
        
        function WU = WaveUp(obj)
            WUPs = zeros(1, length(obj.ConnectedPorts));
            for i=1:length(obj.ConnectedPorts)
                WUPs(i) = WaveUp(obj.ConnectedPorts(i));
            end
            obj.WUPs_utn = WUPs;
            WUPs = [WUPs(1:obj.AdaptedIndex-1), 0, WUPs(obj.AdaptedIndex:end)];
            WU = obj.S(obj.AdaptedIndex,:)*(WUPs');
            obj.WU = WU;
        end
        
        function set.WD(obj, WaveFromParent)
            obj.WD = WaveFromParent;
            WUPs = [obj.WUPs_utn(1:obj.AdaptedIndex-1), WaveFromParent, obj.WUPs_utn(obj.AdaptedIndex:end)];
            b = obj.S*(WUPs');
            for i=1:length(obj.ConnectedPorts)
                if i < obj.AdaptedIndex
                    set(obj.ConnectedPorts(i), 'WD', b(i));
                else
                    set(obj.ConnectedPorts(i), 'WD', b(i+1));
                end
                    
            end
        end
        
    end
    
end

%numNodes;
%numPort;
%numAbsorbedElements;
%R_vect = [Ra, Rb, Rc, Rd, Re, Rf, Rg, Rh, Ri, Rj, Rk];
%R_value = [p1.PortRes, p2.PortRes, p3.PortRes, Rbw.PortRes, Cbw.PortRes, Rout.PortRes, RL.PortRes, C2.PortRes, R2.PortRes, C1.PortRes];

