
startPoint = 1;
resolveAll = 0;

ronscale = 1/10;
cossscale = 1/1e9;

rng(3);
sf = 50;
data = rand(50,2) + [1/sf 0];
data = [data; rand(20,2)./[5 1] + [1/sf 0]];
data = [data; rand(5,2)./2 + .5];
% data(data(:,2) < 50e-12*5e9,:) = [];

data(:,2) = 1./data(:,1)/sf + data(:,2).^3;
data( data(:,2) > 1,:) = [];
data = [data; .95 .95];

figure(1)
hold off;
voronoi(data(:,1), data(:,2));
[vx,vy]=voronoi(data(:,1), data(:,2));
axis square
ylim([0 1])
hold on;
x = 0:.01:1;
plot(x,1./x/sf,':r')

[~,ic] = max(sum(data.^2,2));

hold on;
plot(data(ic,1), data(ic,2),'or')
axis square
xlabel('Ron [pu]')
ylabel('Coss [pu]')

nump = size(data,1);
plabels = arrayfun(@(n) {sprintf('%d', n)}, (1:nump)');
hold on
Hpl = text(data(:,1), data(:,2), plabels, 'FontWeight', ...
      'bold', 'HorizontalAlignment','center', ...
      'BackgroundColor', 'none');


zlims = [70 100];

% plot(vx(:,1),vy(:,1), ':k', 'LineWidth',3)


while(0)
sf = [1 1];
minvec = -data(ic,:).*sf;
for i = 1:20
    %min distance along specified direction to get closer to point k than
    %current point
    ptdist = data-data(ic,:);
    mindist = sum(ptdist.^2,2) ./ sum(ptdist.*minvec,2);

%     mindist = ((data(:,1)-data(ic,1)).^2 + (data(:,2)-data(ic,2)).^2) ./ ...
%         (minvec(1)*(data(:,1)-data(ic,1)) + minvec(2)*(data(:,2)-data(ic,2)));

    mindist(mindist < 0) = inf;
    [~,nextpt] = min(mindist);
    if isinf(mindist(nextpt))
        break
    end
    plot(data(nextpt,1), data(nextpt,2),'d')
    ic = nextpt;
    minvec = -data(ic,:).*sf;
end
end

%%


sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);

modelfile = 'MRBuck'; PLECsModel = 'MRBuck';

        
% find_system(modelfile,'SearchDepth',1, 'IncludeCommented', 'on')
open_system(modelfile,'loadonly');
circuitPath = [modelfile '/' PLECsModel];
set_param(circuitPath,'Commented','on');
clear sim;
simout = sim(modelfile,eps);

for i = 1:length(simout.properties)
    assignin('base',simout.properties{i},eval(['simout.' simout.properties{i}]));
end

set_param(circuitPath,'Commented','off');

Lr = 0;
Lf = 0.25e-6;
% swvec = [0 0; 0 1 ; 0 0; 1 0];
% Ts = .1e-6;
% ts = [.05 .8 .05 .1]*Ts;
swvec = [0 1 ;  1 0];
Ts = [1e-6 .1e-6 5e-6];
Ts = Ts(startPoint);
ts = [ .9  .1]*Ts;
    
% data(:,1) = data(:,1)/10; %ron
% data(:,2) = data(:,2)/1e8; % Coss



sim = SMPSim();
%% Get continuous answer
if(resolveAll)
    nx = 20;
    ny = 21;

    etaSolvedCont = zeros(nx,ny);
    fsSolveCont = zeros(nx,ny);

    %% dfg

    RonVec = linspace(min(data,[],'all'), 1, nx)*ronscale;
    CossVec = linspace(min(data,[],'all'), 1, ny)*cossscale;
    for i = 1:nx
        
        for j = 1:ny
            ron = RonVec(i);
            Coss = CossVec(j);
            CdsL = Coss;
            CdsH = Coss;
            conv = sim.converter;
            top = sim.topology;
            
            top.loadCircuit(circuitPath,swvec,1);
            sim.u = us';
            conv.setSwitchingPattern(1:size(swvec,1), ts)
            
            
            f = @(x) fs_etamin(sim, Vg, Io, x);
            x = fminbnd(f,.01, 10);
        
            newTs = x/1e6;
        
            ts = conv.fullts(1,:)/sum(conv.fullts(1,:))*newTs;
            conv.setSwitchingPattern(1:size(swvec,1), ts)
        
        
%             niter = sim.findValidSteadyState;
            sim.steadyState;
            
            Voloc = find(strcmp(sim.stateNames,'Co'));
            Igloc = find(strcmp(sim.outputNames,'Ig'));
            
            [avgXs, avgYs] = sim.ssAvgs(sim.Xs);
            
            etaSolvedCont(i,j) = avgXs(Voloc)*Io/(Vg*avgYs(Igloc));
            fsSolveCont(i,j) = 1/newTs;

        end
    end

    figure(3)

    etaSolvedValid = etaSolvedCont;
    [RV,CV] = meshgrid(RonVec,CossVec);
    Inval = CV/cossscale < 1./(RV/ronscale)/sf;
    etaSolvedValid = etaSolvedValid.*(~Inval' + 0);

    [xq,yq] = meshgrid(0:.0001:.1, 0:1e-13:2e-10);
    vq = interp2(RonVec,CossVec,etaSolvedCont',xq,yq, 'linear');
    zq = interp2(RonVec,CossVec,log10(fsSolveCont'),xq,yq, 'linear');

    Inval = yq/cossscale < 1./(xq/ronscale)/sf;
    vq = vq.*(~Inval + 0);
    vq(vq==0) = Inf;


%     surf(RonVec,CossVec,etaSolvedValid',fsSolveCont')
    surf(xq*1000,yq*1e12,vq*100,zq, 'LineStyle','none');
    zlim([.8 1])

    %% plot discrete points on top
    figure(3);
    hold on;
    for i = 1:length(data)
        xp = data(i,1)*ronscale;
        yp = data(i,2)*cossscale;
        zp = interp2(RonVec,CossVec,etaSolvedCont',xp,yp, 'linear');
        plot3(xp*1000,yp*1e12,zp*100,'.k','MarkerSize',20)
    end

    box on; grid on;

    xlim([0 100])
    ylim([0 200])
    zlim(zlims)
    xlabel('R{on} [m\Omega]')
    ylabel('C{oss} [pF]')
    zlabel('\eta [%]')

    cb = colorbar;
    clims = [min(log10(fsSolveCont),[],'all') max(log10(fsSolveCont),[],'all')];
    caxis(clims)
    cb.Limits = clims;
    newTicksLin = round(10.^(cb.Ticks)/1e6,1);
    cb.Ticks = log10(newTicksLin*1e6);
    cb.TickLabels = num2str(newTicksLin');
    cb.Label.String = 'fs [MHz]';

    fig  = gcf;
    view(37.5+90,30)
    fig.Units = 'Inches';
    fig.Position = [7.16 .4167 7.16 3.5];
    yl = get(gca,'YLabel');
    set(yl,'rotation',-10,'VerticalAlignment','middle','HorizontalAlignment','right')
    xl = get(gca,'XLabel');
    set(xl,'rotation',20,'VerticalAlignment','middle','HorizontalAlignment','left')
    ytickangle(0)
%     fig.Renderer = 'painters';

end


%% Get ground truth
if(resolveAll)

etaSolved = zeros(length(data),1);
convergeSolve = zeros(length(data),1);
fsSolve = zeros(length(data),1);

for i = 1:length(data)
    ron = data(i,1)*ronscale;
    Coss = data(i,2)*cossscale;
    CdsL = data(i,2)*cossscale;
    CdsH = data(i,2)*cossscale;
    conv = sim.converter;
    top = sim.topology;
    
    top.loadCircuit(circuitPath,swvec,1);
    sim.u = us';
    conv.setSwitchingPattern(1:size(swvec,1), ts)
    
%     sim.steadyState;
    % if(debug)
    %     sim.plotAllStates(1);
    % end
    
    f = @(x) fs_etamin(sim, Vg, Io, x);

%     x = fmincon(f,1,-1,0);
    x = fminbnd(f,.01, 10);

    newTs = x/1e6;

    ts = conv.fullts(1,:)/sum(conv.fullts(1,:))*newTs;
    conv.setSwitchingPattern(1:size(swvec,1), ts)

    sim.steadyState;
%     niter = sim.findValidSteadyState;
    
    Voloc = find(strcmp(sim.stateNames,'Co'));
    Igloc = find(strcmp(sim.outputNames,'Ig'));
    
    [avgXs, avgYs] = sim.ssAvgs(sim.Xs);
    
    etaSolved(i) = avgXs(Voloc)*Io/(Vg*avgYs(Igloc));
%     convergeSolve(i) = niter;
    fsSolve(i) = 1/newTs;
end


end
figure(2)
hold on;


% nBndPts = 25;
% data = [data; 
%     -1*ones(nBndPts,1), linspace(1/nBndPts,1-1/nBndPts,nBndPts)';
%     2*ones(nBndPts,1), linspace(1/nBndPts,1-1/nBndPts,nBndPts)'
%     linspace(1/nBndPts,1-1/nBndPts,nBndPts)', -1*ones(nBndPts,1)
%     linspace(1/nBndPts,1-1/nBndPts,nBndPts)', 2*ones(nBndPts,1)];
origData = data;
data = [data;
    -data(:,1) data(:,2)
    data(:,1) -data(:,2)
    -data(:,1) -data(:,2);
    -data(:,1)+2 data(:,2)
    data(:,1) -data(:,2)+2]



[v,c] = voronoin(data);
data = origData;
% data = data(1:length(data)-4*nBndPts,:);

zdata = etaSolved*100;


cvec =parula;
clength = length(cvec);
cdata = log10(fsSolve);
clims = [min(cdata) max(cdata)];
cdata = min(max(cdata,clims(1)),clims(2));
cdata = cvec(min(256,round((cdata-clims(1))/(diff(clims))*clength)+1),:);
% eColor = (cdata).^2;
eColor = (cdata).^(1/4);

for i = 1:length(data)
    ind = [c{i} c{i}(1)];
    pts = v(ind,:);
    pts = pts(~isinf(sum(pts,2)),:);
%     pts(any(pts<0,2),:) = [];
%     pts(any(pts>1,2),:) = [];
%     if sum(isinf(sum(pts,2))) > 0
%         for j = length(pts):-1:1
%             if isinf(sum(pts(j,:),2))
%                 if j < length(pts)
%                     [ii,jj] = find(vx == pts(j+1,1) & vy ==pts(j+1,2));
%                     ii = 2-ii+1;
%                     inds = sub2ind(size(vx),ii,jj);
%                     newpts = [vx(inds) vy(inds)];
%                 else
%                     newpts = [];
%                 end
%                 if j > 1
%                     [ii,jj] = find(vx == pts(j-1,1) & vy ==pts(j-1,2));
%                     ii = 2-ii+1;
%                     inds = sub2ind(size(vx),ii,jj);
%                     newpts = [newpts; vx(inds) vy(inds)];
%                 end
%                 pts = [pts(1:max(1,j-1),:); newpts; pts(min(length(pts),j+1),:)];
%             end
%         end
%     end


    pts(:,1) = pts(:,1)*ronscale*1000;
    pts(:,2) = pts(:,2)*cossscale*1e12;
    if numel(pts) < 6
        continue
    end
    patch(pts(:,1), pts(:,2), zdata(i)*ones(1,length(pts)),cdata(i,:),'EdgeColor','none');
%     patch(pts(:,1), pts(:,2), zdata(i)*ones(1,length(pts)),cdata(i,:),'EdgeColor',eColor(i,:));
    for j = 1:length(pts)-1
        ind2 = [ j j+1 j+1 j];
%         patch(v(ind(ind2),1), v(ind(ind2),2), [zdata(i)*ones(1,2) zeros(1,2)],cdata(i,:),'EdgeColor',(cdata(i,:)).^2);
        patch(pts(ind2,1), pts(ind2,2), [zdata(i)*ones(1,2) zeros(1,2)],cdata(i,:),'EdgeColor','none');
%         patch(pts(ind2,1), pts(ind2,2), [zdata(i)*ones(1,2) zeros(1,2)],cdata(i,:),'EdgeColor',eColor(i,:));

    end
end
xlim([0 100])
ylim([0 1000])
zlim(zlims)
xlabel('R{on} [m\Omega]')
ylabel('C{oss} [pF]')
zlabel('\eta [%]')
box on; grid on;

cb = colorbar('NorthOutside');
caxis(clims)
cb.Limits = clims;
newTicksLin = round(10.^cb.Ticks/1e6,1);
cb.Ticks = log10(newTicksLin*1e6);
cb.TickLabels = num2str(newTicksLin');
cb.Label.String = 'fs [MHz]';

alpha(0.3)



figure(1)
% [~,ic] = max(mean(data,2).^2);
ic = [68 19 5];
ic = ic(startPoint)
oldpt = [];

% allPts = zeros(100,2);
% 
% allPts(1,:) = data(ic,:);

%% Analyze circuit
for i = 1:100
    ron = data(ic,1)*ronscale;
    Coss = data(ic,2)*cossscale;
    CdsL = data(ic,2)*cossscale;
    CdsH = data(ic,2)*cossscale;
    sim = SMPSim();
    conv = sim.converter;
    top = sim.topology;
   
    top.loadCircuit(circuitPath,swvec,1);
    sim.u = us';
    conv.setSwitchingPattern(1:size(swvec,1), ts)
    
    sim.steadyState;
    % if(debug)
    %     sim.plotAllStates(1);
    % end
    
%     niter = sim.findValidSteadyState;
    
%     Igloc = find(strcmp(sim.stateNames,'L1'));
    Voloc = find(strcmp(sim.stateNames,'Co'));
    Igloc = find(strcmp(sim.outputNames,'Ig'));
    
    [avgXs, avgYs] = sim.ssAvgs(sim.Xs);
    
    eta = avgXs(Voloc)*Io/(Vg*avgYs(Igloc));

    %% perturb ron
        ron = ron*1.05;
        top.loadCircuit(circuitPath,swvec,1);
        sim.u = us';
        conv.setSwitchingPattern(1:size(swvec,1), ts)
        sim.steadyState;
%         niter = sim.findValidSteadyState;
        [avgXs, avgYs] = sim.ssAvgs(sim.Xs);
        
        etaRon = avgXs(Voloc)*Io/(Vg*avgYs(Igloc));
    
    %% perturb Coss
        ron = ron/1.05;
        CdsH = Coss*1.2;
        CdsL = CdsH;
        top.loadCircuit(circuitPath,swvec,1);
        sim.u = us';
        conv.setSwitchingPattern(1:size(swvec,1), ts)
        sim.steadyState;
%         niter = sim.findValidSteadyState;
        [avgXs, avgYs] = sim.ssAvgs(sim.Xs);
        
        etaCoss = avgXs(Voloc)*Io/(Vg*avgYs(Igloc));

        minvec = [(etaRon-eta)/(data(ic,1)*.05) (etaCoss-eta)/(data(ic,2)*.2) ];

   %% perturb fs
        CdsH = Coss;
        CdsL = Coss;
        top.loadCircuit(circuitPath,swvec,1);
        sim.u = us';
        conv.setSwitchingPattern(1:size(swvec,1), ts*1.05)
        sim.steadyState;
%         niter = sim.findValidSteadyState;
        [avgXs, avgYs] = sim.ssAvgs(sim.Xs);
        
        etaTs = avgXs(Voloc)*Io/(Vg*avgYs(Igloc));

    minvec = [(etaRon-eta)/(data(ic,1)*.05) (etaCoss-eta)/(data(ic,2)*.2) (etaTs-eta)/(.05)];

    
    
    
    
    ptdist = data-data(ic,:);
    mindist = sum(ptdist.^2,2) ./ sum(ptdist.*minvec(1:2),2);
    mindist(mindist < 0) = inf;
    [deltaDist,nextpt] = mink(mindist,2);
    if nextpt(1) == oldpt
        nextpt = nextpt(2);
        deltaDist = deltaDist(2);
    else
        nextpt = nextpt(1);
        deltaDist = deltaDist(1);
    end
    if isinf(mindist(nextpt))
        dT = (1+.05*minvec(3)*10);
        dT = min(max(dT,.5),1.5);
        ts = ts*dT;
        continue
    end

    dT = (1+deltaDist*minvec(3)*10);
    dT = min(max(dT,.5),1.5);
    ts = ts*dT;



    p1 = data(ic,:);
    p2 = data(ic,:) + deltaDist*minvec(1:2)/2;%minvec(1:2)/norm(minvec(1:2),2)/10;
    dp = p2-p1;
    q = quiver(p1(1),p1(2),dp(1),dp(2),0,'Linewidth',2);
    q.MaxHeadSize = 0.5;
%     annotation('textarrow',[data(ic,1) data(ic,1)+minvec(1)/norm(minvec,2)/10],[data(ic,2) data(ic,2)+minvec(2)/norm(minvec,2)/10],'String','test ','FontSize',13,'Linewidth',2)
    plot(data(nextpt,1), data(nextpt,2),'d','LineWidth',3)
    [nextpt eta sum(ts)*1e6]
    oldpt = ic;
    ic = nextpt;


    
    if ~(nextpt == oldpt)
        figure(2)
        plot3(data(oldpt,1)*ronscale*1000, data(oldpt,2)*cossscale*1e12, etaSolved(oldpt)*100, '.k','LineWidth',3,'MarkerSize',10)
        q = quiver3(data(oldpt,1)*ronscale*1000, data(oldpt,2)*cossscale*1e12, etaSolved(oldpt)*100,dp(1)*ronscale*1000,dp(2)*cossscale*1e12,0,'k','Linewidth',2);
        q.MaxHeadSize = 0;%3/sqrt(sum(dp.^2));

%         p3 = p1 + deltaDist/2.*minvec(1:2);
    
        traj = [ p2; p2; data(nextpt,:)];
        height = [etaSolved(oldpt) etaSolved(ic) etaSolved(ic)];
        plot3(traj(:,1)*ronscale*1000, traj(:,2)*cossscale*1e12, height*100, ':r','LineWidth',2)
    
        figure(1)
    end


end

figure(2)

plot3(data(ic,1)*ronscale*1000, data(ic,2)*cossscale*1e12, etaSolved(ic)*100, '.r','LineWidth',3,'MarkerSize',10)
fig  = gcf;
view(37.5+90,30)
fig.Units = 'Inches';
fig.Position = [0 .4167 3.5 3.5];
yl = get(gca,'YLabel');
set(yl,'rotation',-10,'VerticalAlignment','middle','HorizontalAlignment','right')
xl = get(gca,'XLabel');
set(xl,'rotation',20,'VerticalAlignment','middle','HorizontalAlignment','left')
ytickangle(0)
% fig.Renderer = 'painters';


function sigma = fs_etamin (sim, Vg, Io, Ts1u)
    Ts = Ts1u/1e6; 
    conv = sim.converter;
    swvec = conv.swvec;
    ts = conv.fullts(1,:)/sum(conv.fullts(1,:))*Ts;
    conv.setSwitchingPattern(1:size(swvec,1), ts)
    
    sim.steadyState;
    % if(debug)
    %     sim.plotAllStates(1);
    % end
    
%     niter = sim.findValidSteadyState;
    
%     Igloc = find(strcmp(sim.stateNames,'L1'));
    Voloc = strcmp(sim.stateNames,'Co');
    Igloc = strcmp(sim.outputNames,'Ig');
    
    [avgXs, avgYs] = sim.ssAvgs(sim.Xs);
    
    eta = avgXs(Voloc)*Io/(Vg*avgYs(Igloc));
    sigma = 1-eta;

    if sigma < 0 
        sigma = .5;
    end
end
