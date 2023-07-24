classdef MOSFET < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ron
        qg
        Rjc
        CossVds
        IsdVsd
        RonT
        Vds
        CeqQ
        CeqE
        Vbr
    end
    
    methods
        function obj = MOSFET(ron, qg, Rjc, CossVds, IsdVsd, RonT, Vbr)
            obj.ron = ron;
            obj.qg = qg;
            obj.Rjc = Rjc;
            obj.CossVds = CossVds;
            obj.IsdVsd = IsdVsd;
            obj.RonT = [RonT(:,1) RonT(:,2)*ron];
            obj.Vbr = Vbr;
            
            Vds = linspace(0,Vbr,1000);

            Coss = interp1(CossVds(:,1), CossVds(:,2), Vds);
            CeqE = 2./Vds.^2.*cumtrapz(Vds, Coss.*Vds);
            CeqE(1) = Coss(1);
            CeqQ = 1./Vds.*cumtrapz(Vds, Coss);
            CeqQ(1) = Coss(1);
            obj.Vds = Vds;
            obj.CeqQ = CeqQ;
            obj.CeqE = CeqE;
        end
        
        function updateTemp(obj, T)
            obj.ron = min(interp1(obj.RonT(:,1), obj.RonT(:,2), T, 'linear', 'extrap'),obj.RonT(end,2));
        end
    end
    
end

