classdef Diode < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Cr
        Vr
        Vf
        Rd
        Vbr
        Ifmax
        T
        CeqQ
        CeqE
        VfIdT
        Rja
    end
    
    methods
        function obj = Diode(Vbr, CrVr, VfIdT, Rd, Rja)
            obj.Vbr = Vbr;
            obj.Rja = Rja;
            
            if(numel(CrVr) == 1)
                obj.Cr = CrVr;
                obj.Vr = 0;
            elseif(min(size(CrVr)) == 2)
%                 dim = find(size(Cr) == 2);
                Vr = linspace(0,Vbr,1000);
                Cd = interp1(CrVr(:,1), CrVr(:,2), Vr, 'linear', 'extrap');
                CeqE = 2./Vr.^2.*cumtrapz(Vr, Cd.*Vr);
                CeqE(1) = Cd(1);
                CeqQ = 1./Vr.*cumtrapz(Vr, Cd);
                CeqQ(1) = Cd(1);
                obj.Vr = Vr;
                obj.CeqQ = CeqQ;
                obj.CeqE = CeqE;
            end
            
            if(numel(VfIdT) == 1)
                obj.Vf = Vf;
                obj.Rd = Rd;
            else
%                 dim = find(size(Cr) == 2);
                Vf_Id_T = scatteredInterpolant(VfIdT(:,2), VfIdT(:,3), VfIdT(:,1));
                obj.VfIdT = Vf_Id_T;
                obj.Vf = 0;
                obj.Rd = 0;
            end
        end
        
        function Vf = diodeVoltage(obj, Id, T)
            if(obj.Vf == 0  && obj.Rd==0)
                if(nargin == 3)
                    Vf = obj.VfIdT(Id, T); 
                else
                    Vf = obj.VfIdT(Id, 25*ones(size(Id))); 
                end
            else
                Vf = obj.Vf + Id*obj.Rd;
            end
        end
    end
    
end

