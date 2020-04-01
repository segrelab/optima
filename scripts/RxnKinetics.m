classdef RxnKinetics
    %RXNKINETICS A class to contain reaction equations and ODE functions
    %for extracellular reactions
    
    properties
        %Property1
    end
    
    methods
        function obj = RxnKinetics()
        end
        
        %         function obj = RxnKinetics(inputArg1,inputArg2)
        %             %RXNKINETICS Construct an instance of this class
        %             %   Detailed explanation goes here
        %             obj.Property1 = inputArg1 + inputArg2;
        %         end
        %
        %         function outputArg = method1(obj,inputArg)
        %             %METHOD1 Summary of this method goes here
        %             %   Detailed explanation goes here
        %             outputArg = obj.Property1 + inputArg;
        %         end
    end
    methods (Static)
        %Enzymatic reaction
        %Args: enzyme concentration, substrate concentration,
        %catalysis rate (s^-1), Michaelis constant
        function res = enz(e,s,kcat,km)
            res = RxnKinetics.mm(e,s,kcat,km);
        end
        
        %Enzymatic reaction with product inhibition
        %Args: enzyme concentration, substrate concentration, product
        %concentration, catalysis rate (s^-1), Michaelis constant,
        %inhibition constant
        function res = enzpi(e,s,p,kcat,km,kp)
            res = (kcat * e * s) / ((km * (1 + (p/kp))) + s);
        end
        
        %basic michaelis-menten sans inhibition
        function res = mm(e,s,kcat,km)
            res = (kcat * e * s)/ (km + s);
        end
        
        %ODE timestep functions for a enzymatic reaction with product
        %inhibition
        %if the final arg is empty or 0, it will not have product
        %inhibition. If given, the final argument is the inhibition
        %constant
        function res = enzdydt(e,s,p,kcat,km,varargin)
            res = RxnKinetics.mmdydt(e,s,p,kcat,km,varargin);
        end
        function res = mmdydt(e,s,p,kcat,km,varargin)
            ki = 0;
            if (~isempty(varargin))
                ki = varargin{1};
            end
            
            if (ki == 0)
                ds = RxnKinetics.mm(e,s,kcat,km);
            else
                ds = RxnKinetics.enzpi(e,s,p,kcat,km,ki);
            end
            %don't let substrate drain into the negative
            if (ds > s)
                ds = s;
            end
            res = [0;
                -ds;
                ds];
        end
        
        %ODE function for describing an enzymatic reaction as a series with
        %intermediate species, no explicit inhibition. Includes reversible 
        %binding reactions
        %E + S <-> ES <-> EP <-> E + P
        %make r2 irreversible by setting k2r=0
        function res = dydtintermediate(e,s,p,es,ep,k1f,k1r,k2f,k2r,k3f,k3r)
            res = [-(e * s * k1f) + (es * k1r) + (ep * k3f) - (e * p * k3r); %e
                -(e * s * k1f) + (es * k1r); %s
                (ep * k3f) - (e * p * k3r); %p
                (e * s * k1f) - (es * k1r) - (es * k2f) + (ep * k2r); %es
                (es * k2f) - (ep * k2r) - (ep * k3f) + (e * p * k3r)]; %ep
        end
        
        %ODE function for modeling an enzymatic reaction as two reactions
        %catalyzed by different enzymes, with product inhibition
        %E + S --e1> EP --e2> E + P
        %set either Ki to 0 to disable product inhibition
        function res = dydt2step(e1,e2,s1,s2,p,kcat1,km1,ki1,kcat2,km2,ki2)
            res = [0;
                0;
                -RxnKinetics.enzpi(e1,s1,s2,kcat1,km1,ki1);
                RxnKinetics.enzpi(e1,s1,s2,kcat1,km1,ki1) - RxnKinetics.enzpi(e2,s2,p,kcat2,km2,ki2);
                RxnKinetics.enzpi(e2,s2,p,kcat2,km2,ki2)];
        end
        
        
    end
end

