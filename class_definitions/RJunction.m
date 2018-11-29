classdef RJunction < WDF
    %RJUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        S;
        WU;
        WD;
        Ports;
    end
    
    methods
        function obj = RJunction(MnaMatrix, R_vect, R_value, adaptedElement)
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
            
            %compute scattering matrix: S = I + 2[0 R] * X^(-1) * [0 I]'
            R_diag = diag(R_vect);
            left_term = 2*[zeros(numPort, (numNodes-1) ), R_diag, zeros(numPort, numAbsorbedElements)];
            right_term = [zeros(numPort, (numNodes-1) ), eye(numPort), zeros(numPort, numAbsorbedElements)];
            S = eye(numPort) + left_term*(X_inv*(right_term'));
            
            [R_PortRes, param, cond] = solve(S(1,1)==0, Ra, 'ReturnConditions', true);
            assume(cond);
            R_PortRes = double(subs(R_PortRes, R_vect(2:end), R_value));
            R_value = [R_PortRes, R_value];
            obj.S = double(subs(S, R_vect, R_value));
            
        end
        
    end
    
end

%numNodes;
%numPort;
%numAbsorbedElements;
%R_vect = [Ra, Rb, Rc, Rd, Re, Rf, Rg, Rh, Ri, Rj, Rk];
%R_value = [p1.PortRes, p2.PortRes, p3.PortRes, Rbw.PortRes, Cbw.PortRes, Rout.PortRes, RL.PortRes, C2.PortRes, R2.PortRes, C1.PortRes];

