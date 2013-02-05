classdef CKalmanUdHess < handle
    %CKALMAN 
    
    properties
        %Scalar
        nx
        T
        
        %Vectors
        Xest
        Xextr
        Xbase
        Dksi
        
        % Matrix
        F
        G
        
        Dest
        Dextr
    end
    
    methods
        function K = CKalmanUdHess(X0, Dest0, F, G, Dksi0)
            K.Xest = X0;
            K.F = F; 
            K.Dest = Dest0;
            K.G = G;
            K.Dksi = Dksi0;
            K.nx = length(X0);
        end
        
        function Xextrap = Extrapolate(K)
            K.Xextr = K.F*K.Xest;
            K.Dextr = K.F* K.Dest * (K.F)' + K.G*K.Dksi*K.G';
            Xextrap = K.Xextr;
        end
        
        function setBase(K, Xb)
            K.Xbase = Xb;
        end
        
        function Estimate(K, u, Hess)
            K.Dest = inv( inv(K.Dextr) - Hess );
            K.Xest = K.Xextr + K.Dest*(u + Hess*(K.Xextr - K.Xbase));
        end
        
        function xest = getXest(K)
            xest = K.Xest;
        end
    end
end
