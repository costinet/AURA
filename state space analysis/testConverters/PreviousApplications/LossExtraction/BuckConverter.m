classdef BuckConverter < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        V
        C
        L
        Rl
        Rl_400k
        Rl_800k
        fs
        Ts
        dt
        Vdr
        HS_FET
        LS_FET
        LS_Diode
        
        As
        Bs
        
        Cw
    end
    
    methods
        function obj = BuckConverter(V, C, L, Rl, Rl_400k, Rl_800k, fs, Ts, dt, Vdr, HS_FET, LS_FET, As, Bs, Cw, LS_Diode)
            obj.V = V;
            obj.C = C;
            obj.L = L;
            obj.Rl = Rl;
            obj.Rl_400k = Rl_400k;
            obj.Rl_800k = Rl_800k;
            obj.fs = fs;
            obj.Ts = Ts;
            obj.dt = dt;
            obj.Vdr = Vdr;
            obj.HS_FET = HS_FET;
            obj.LS_FET = LS_FET;

            obj.As = As;
            obj.Bs = Bs;
            if(nargin > 14)
                obj.Cw = Cw;
            end
            if(nargin >15)
                obj.LS_Diode = LS_Diode;
            end
        end
        
        
        
        function [Xss, ys, t, ts] = simulate(obj, Io, Vg)
            Dapprox = obj.V./Vg;
            u = [Vg Io]';
            ts = [Dapprox*obj.Ts - obj.dt,obj.dt,(1-Dapprox)*obj.Ts - obj.dt, obj.dt];
%             ts = [Dapprox*obj.Ts - obj.dt,obj.dt,(1-Dapprox)*obj.Ts - obj.dt, 5.1764e-09];
            if sum(ts<0)
                errLoc = find(ts<0,1,'first');
                ts(mod(errLoc+2,4)) = ts(mod(errLoc+2,4)) + ts(errLoc);
                ts(errLoc) = 0;
            end
                
            
            Cp = (obj.HS_FET.CeqQ(find(obj.HS_FET.Vds>=Vg,1,'first'))+ obj.HS_FET.CeqQ(find(obj.HS_FET.Vds>=Vg,1,'first')))*1e-12;
            
            if(obj.Cw > 0)
                K = [Cp 0 0 0; 0 obj.L 0 0; 0 0 obj.C 0; 0 0 0 obj.Cw];
            else
                K = [Cp 0 0 ; 0 obj.L 0 ; 0 0 obj.C];
            end

            for i = 1:size(obj.As,3)
                As(:,:,i) = K^-1*obj.As(:,:,i);
                Bs(:,:,i) = K^-1*obj.Bs(:,:,i);
            end

            [ Xss] = SS_Soln( As, Bs, ts, u);

            %% Model Diode nonlinearities and Output Regulation

            diodcon2 = (Xss(1,3) < 0);
            diodcon4 = (Xss(1,5) > Vg) || (Xss(1,5) < -2);
            Voerror = abs(obj.V- sum(Xss(3,1:end-1).*ts)/obj.Ts) > 1;
            
            errVec = [diodcon2 diodcon4 Voerror]

            while(sum(errVec));
                if diodcon2
                    dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 2, 1e-9, 3);
                    dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 2, -1e-9, 3);

                    dXds = (dXsp - dXsn)/2e-9;
                    deltaT = max(min(Xss(1,3)/dXds(1,3), ts(2)),0);
                    ts(2) = ts(2) - deltaT;
                    ts(3) = ts(3) + deltaT;
                end

                [ Xss] = SS_Soln( As, Bs, ts, u);

                if diodcon4
                    dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 4, 1e-9, 1);
                    dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 4, -1e-9, 1);

                    dXds = (dXsp - dXsn)/2e-9;
                    deltaT = max(min((Xss(1,5)-Vg)/dXds(1,5), ts(4)),0);
                    ts(4) = ts(4) - deltaT;
                    ts(1) = ts(1) + deltaT;
                end



                [ Xss] = SS_Soln( As, Bs, ts, u);

                if Voerror
                    dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 3, 1e-9, 1);
                    dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 3, -1e-9, 1);

                    dXds = (dXsp - dXsn)/2e-9;
                    Voavg = sum(Xss(3,1:end-1).*ts)/obj.Ts;
                    dVdt = mean(dXds(3,:));
                    deltaT = max(min((Voavg-obj.V)/dVdt, ts(3)),-ts(1));
                    ts(3) = ts(3) - deltaT;
                    ts(1) = ts(1) + deltaT;
                end

                [ Xss] = SS_Soln( As, Bs, ts, u);
                
                

                diodcon2 = (Xss(1,3) < -.5);
                diodcon4 = (Xss(1,5) > Vg+.5) || (Xss(1,5) < -.5);
                Voerror = abs(obj.V- sum(Xss(3,1:end-1).*ts)/obj.Ts) > 1;
                
                errVec = [diodcon2 diodcon4 Voerror]
            end

            [ ys, t ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u );
            
%             figure(50)
%             plot(t, ys(1,:));
%             hold on;
%             plot(sum(ts(1:2)), Xss(1,3),'o');
%             plot(sum(ts(1:4)), Xss(1,5),'o');
%             qwew=1;
        end
        
        
        function [ Pcond, Pg, Pbd, Poss, Pq, Pov, Pboot, Pcore, Voavg, Temps] = CalculateLosses(obj, ys, t, ts, Xss, Vg, Io) 
            iL = ys(2,:);
            Vp = ys(1,:);
            Vo = ys(3,:);

            Voavg = mean(Vo);


            %% Loss Model
            rgp_HS = 2.5+.95;
            rgn_HS = 1.5+1.5;

            irms = sqrt(mean(iL.^2));
            irms1 = sqrt(mean(iL(Vp>Vg-1).^2));
            irms2 = sqrt(mean(iL(Vp<1).^2));
            iacrms = sqrt(mean((iL-mean(iL)).^2));

            %% Constant Losses
            Pq = 50e-6*obj.V + Vg^2/490e3+.1;
            Pboot = ones(1,length(Io)).*(3e-12*.5*Vg^2*obj.fs + 1*obj.Vdr*obj.HS_FET.qg*obj.fs/obj.Vdr);
            Pg = ones(1,length(Io)).*(obj.Vdr*obj.LS_FET.qg + obj.Vdr*obj.HS_FET.qg)*obj.fs;
            Pcore = 0.342;

            %% Coss Losses
            CeqE = (obj.HS_FET.CeqE(obj.HS_FET.Vds<=Vg) + fliplr(obj.LS_FET.CeqE(obj.HS_FET.Vds<=Vg)))*1e-12 + ...
                obj.LS_Diode.CeqE(obj.LS_Diode.Vr<=Vg) + 10e-12;
            Vp3 = Xss(1,3);
            Poss1 = .5*CeqE(min(find(obj.HS_FET.Vds>Vp3,1),length(CeqE)))*Vp3^2*obj.fs;
            Vp5 = Xss(1,5);
            Poss2 = .5*CeqE(end)*Vg^2*obj.fs - .5*CeqE(min(find(obj.HS_FET.Vds>Vp5,1),length(CeqE)))*Vp5^2*obj.fs;
            Poss = Poss1 + Poss2;

            %% Overlap Losses
            Imin = min(iL);
            Imax = max(iL);
            Pov = .5*(1.4*rgp_HS*obj.HS_FET.qg/obj.Vdr)*Imin*Vg*obj.fs + ...
                .5*(1.4*rgn_HS*obj.HS_FET.qg/obj.Vdr)*Imax*Vg*obj.fs ;

            %% Temperature-dependent
            T_LsD = 25;
            T_HsM = 25;
            T_LsM = 25;

            dTs = [10 10 10]';

            while(sum(abs(dTs)) > 3)

                %% Conduction Losses
                Pcond_LS = irms1.^2.*obj.LS_FET.ron;
                Pcond_HS = irms2.^2.*obj.HS_FET.ron;
                Pcond_L = irms.^2*obj.Rl + iacrms^2*obj.Rl_400k;
                Pcond = Pcond_LS + Pcond_HS + Pcond_L;

                %% Dead Time Body Diode Conduction
                tbd_ls = find(t<ts(1)+obj.dt & t>ts(1)+ts(2));
                tbd_hs = find(t<sum(ts(1:3))+obj.dt & t>sum(ts(1:3)));
                Vsd_HS = interp1(obj.HS_FET.IsdVsd(:,2), obj.HS_FET.IsdVsd(:,1), max(-iL(tbd_hs),0));
                Vsd_LS = diodeVoltage(obj.LS_Diode, max(iL(tbd_ls),0), T_LsD*ones(size(tbd_ls)));
                Ebd_LS = 0;
                Ebd_HS = 0;
                if(length(tbd_ls)>1 )
                    Ebd_LS = trapz(t(tbd_ls), iL(tbd_ls).*Vsd_LS);
                end
                if(length(tbd_hs)>1)
                    Ebd_HS = trapz(t(tbd_hs), iL(tbd_hs).*Vsd_HS );
                end
                Pbd = (Ebd_LS+Ebd_HS)*obj.fs;

                P_LsM = Pcond_LS + Poss1;
                P_HsM = Pcond_HS + Ebd_HS*obj.fs + Poss2 + Pov;
                P_LsD = Ebd_LS*obj.fs;

                dTs = [T_LsD- (P_LsD*obj.LS_Diode.Rja + 27);
                    T_LsM - (P_LsM*obj.LS_FET.Rjc + 27);
                    T_HsM - (P_HsM*obj.HS_FET.Rjc + 27)];

                T_HsM = P_HsM*obj.HS_FET.Rjc + 27;
                T_LsM = P_LsM*obj.LS_FET.Rjc + 27;
                T_LsD = P_LsD*obj.LS_Diode.Rja + 27;

                Temps = [T_HsM T_LsM T_LsD];

                updateTemp(obj.LS_FET, T_LsM);
                updateTemp(obj.HS_FET, T_HsM);

%                 [T_HsM obj.HS_FET.ron]
            end
        end
        
         function [ Pcond, Pconst, Poss, Esw, Ion, Ioff] = FindEsw(obj, ys, t, ts, Xss, Vg, Io, cs, ks) 
            iL = ys(2,:);
            Vp = ys(1,:);
            Vo = ys(3,:);

            Voavg = mean(Vo);


            %% Loss Model
            irms = sqrt(mean(iL.^2));
            irms1 = sqrt(mean(iL(Vp>Vg-1).^2));
            irms2 = sqrt(mean(iL(Vp<1).^2));
            iacrms = sqrt(mean((iL-mean(iL)).^2));

            %% Constant Losses
            Pq = 50e-6*obj.V + Vg^2/490e3+.1;
            Pboot = ones(1,length(Io)).*(3e-12*.5*Vg^2*obj.fs + 1*obj.Vdr*obj.HS_FET.qg*obj.fs/obj.Vdr);
            Pg = ones(1,length(Io)).*(obj.Vdr*obj.LS_FET.qg + obj.Vdr*obj.HS_FET.qg)*obj.fs;
            Pcore = 0.342;  %.342
            Pconst = Pcore + Pg + Pboot + Pq;

            %% Coss Losses
            Cpara = 5e-12;
            CeqE = (obj.HS_FET.CeqE(obj.HS_FET.Vds<=Vg) + fliplr(obj.LS_FET.CeqE(obj.HS_FET.Vds<=Vg)))*1e-12 + ...
                obj.LS_Diode.CeqE(obj.LS_Diode.Vr<=Vg) + Cpara;
            Vp3 = Xss(1,3);
            Poss1 = .5*CeqE(min(find(obj.HS_FET.Vds>Vp3,1),length(CeqE)))*Vp3^2*obj.fs;
            Vp5 = Xss(1,5);
            Poss2 = .5*CeqE(end)*Vg^2*obj.fs - .5*CeqE(min(find(obj.HS_FET.Vds>Vp5,1),length(CeqE)))*Vp5^2*obj.fs;
            Poss = Poss1 + Poss2;

            %% Temperature-dependent
            T_LsD = 25;
            T_HsM = 25;
            T_LsM = 25;

            dTs = [10 10 10]';

            while(sum(abs(dTs)) > 3)

                %% Conduction Losses
                Pcond_LS = irms1.^2.*obj.LS_FET.ron;
                Pcond_HS = irms2.^2.*obj.HS_FET.ron;
                Pcond_L = irms.^2*obj.Rl + iacrms^2*obj.Rl_400k;
                Pcond = Pcond_LS + Pcond_HS + Pcond_L;

                %% Dead Time Body Diode Conduction
                tbd_ls = find(t<ts(1)+obj.dt & t>ts(1)+ts(2));
                tbd_hs = find(t<sum(ts(1:3))+obj.dt & t>sum(ts(1:3)));
                Vsd_HS = interp1(obj.HS_FET.IsdVsd(:,2), obj.HS_FET.IsdVsd(:,1), max(-iL(tbd_hs),0));
                Vsd_LS = diodeVoltage(obj.LS_Diode, max(iL(tbd_ls),0), T_LsD*ones(size(tbd_ls)));
                Ebd_LS = 0;
                Ebd_HS = 0;
                if(length(tbd_ls)>1 )
                    Ebd_LS = trapz(t(tbd_ls), iL(tbd_ls).*Vsd_LS);
                end
                if(length(tbd_hs)>1)
                    Ebd_HS = trapz(t(tbd_hs), iL(tbd_hs).*Vsd_HS );
                end
                Pbd = (Ebd_LS+Ebd_HS)*obj.fs;
                
                %% Esw
                
                Ion = Xss(2,1);
                Ioff = Xss(2,2);
                LS_Ioff =  Xss(2,4);
                
                Vdc = Vg;
                [HSEon, HSEoff] = solveEsw(obj, Ion, Ioff, LS_Ioff, cs, ks, Vdc);
                Esw = HSEon + HSEoff;

                P_LsM = Pcond_LS + Poss1;
                P_HsM = Pcond_HS + Ebd_HS*obj.fs + Esw*obj.fs*1e-6;
                P_LsD = Ebd_LS*obj.fs;

                dTs = [T_LsD- (P_LsD*obj.LS_Diode.Rja + 27);
                    T_LsM - (P_LsM*obj.LS_FET.Rjc + 27);
                    T_HsM - (P_HsM*obj.HS_FET.Rjc + 27)];

                T_HsM = P_HsM*obj.HS_FET.Rjc + 27;
                T_LsM = P_LsM*obj.LS_FET.Rjc + 27;
                T_LsD = P_LsD*obj.LS_Diode.Rja + 27;

                Temps = [T_HsM T_LsM T_LsD];

                updateTemp(obj.LS_FET, T_LsM);
                updateTemp(obj.HS_FET, T_HsM);


            end
         end
        
         function [Psw, Esw] = fitDPT(obj, ys, param, prewarp)
            Ion = ys(2,:);
            Ioff= ys(1,:);
            Vg = ys(3,:);
            
            [cs, ks] = getCKs(obj, param, prewarp);
            
            Vdc = Vg;
            LS_Ioff = Ion;

            [HSEon, HSEoff] = solveEsw(obj, Ion, Ioff, LS_Ioff, cs, ks, Vdc);
            Esw = HSEon + HSEoff;
            
            Psw = Esw*obj.fs*1e-6;
         end
         
         
         function [HSEon, HSEoff] = solveEsw(obj, Ion, Ioff, LS_Ioff, cs, ks, Vdc)
             %% Esw
            
            HSEon = (cs(1) + cs(2)*Ion + cs(3)*Ion.^2).*(ks(1)*Vdc + ks(2)*Vdc.^2) + ks(5);
            HSEoff = (cs(4) + cs(5)*Ioff + cs(6)*Ioff.^2).*(ks(3)*Vdc + ks(4)*Vdc.^2);
%             HSEoff = (cs(4) + cs(5)*Ioff).*(ks(3)*Vdc + ks(4)*Vdc.^2);
              
            ssind = find(Ion<0 & LS_Ioff <0);
            HSEon(ssind) = (cs(7)*abs(Ion(ssind)) + cs(8)*Ion(ssind).^2).*(ks(6)*Vdc(ssind) + ks(7)*Vdc(ssind).^2);
            
            not_valid = find(Ioff<3.5);
            HSEoff(not_valid) = 0;
            %(cs(4) + cs(5)*abs(LS_Ioff(ssind)) + cs(6)*LS_Ioff(ssind).^2).*(ks(3)*Vdc(ssind) + ks(4)*Vdc(ssind).^2);%+ ks(5)/3;  %% Actually LS Eoff
         
                      
         end
         
         function [cs, ks] = getCKs(obj, param, prewarp)
            param = param.*prewarp;

            cs = param(1:8);
            ks = param(9:15); 
         end
    end
    
end

