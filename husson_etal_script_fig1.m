
% matlab script for plotting Fig. 1 in Husson et al. (2022)

load heat_content_u100m.mat
load indgeo
load Sea_ice_data


% June SIE using AO sea ice concentration data
ICEconc_meanJune19791999=mean(ICEconc(:,:,6,find(x_ice==1979):find(x_ice==1999)),4);
ICEconc_meanJune19992009=mean(ICEconc(:,:,6,find(x_ice==1999):find(x_ice==2009)),4);
ICEconc_meanJune20092017=mean(ICEconc(:,:,6,find(x_ice==2009):find(x_ice==2017)),4);
ICEconc_meanJune2006=mean(ICEconc(:,:,6,x_ice==2006),4);
ICEconc_meanJune2012=mean(ICEconc(:,:,6,x_ice==2012),4);
ICEconc_meanJune2016=mean(ICEconc(:,:,6,x_ice==2016),4);


x=year;

% Define the reference period
refPeriodString='1999to2009';
meanPeriod_start_year=1999; % definerer når referanseperioden skal starte
meanPeriod_end_year=2009; % definerer når referanseperioden skal slutte
x_start=find(x==meanPeriod_start_year); 
x_end=find(x==meanPeriod_end_year);
ant_aar=x(x_end)-x(x_start);

% Heat content mean
Q_mean_refPeriod=nan(84,80);
for i=1:I
    for j=1:J
        tmp_Q=squeeze(Q_u100m(i,j,x_start:x_end));
        ind_fin=isfinite(tmp_Q);
         if sum(ind_fin)>ant_aar*0.66 && sum(ind_fin(1:5))>1 && sum(ind_fin(6:end))>1 % Demand for data amount and spread over the time period
            Q_mean_refPeriod(i,j)=mean(tmp_Q(ind_fin))*10^(-6);
         end
    end
end
clear i j p S ind_fin

% Heat content ano
Q_2006=squeeze(Q_u100m(:,:,x==2006))*10^(-6);
Q_2012=squeeze(Q_u100m(:,:,x==2012))*10^(-6);
Q_2016=squeeze(Q_u100m(:,:,x==2016))*10^(-6);
Q_ano2006=Q_2006-Q_mean_refPeriod;
Q_ano2012=Q_2012-Q_mean_refPeriod;
Q_ano2016=Q_2016-Q_mean_refPeriod;

% plot mean
figure, hold on, m_proj('Stereographic','lon',38,'lat',76,'rad',[18 68],'rec','on'); % Hele Barentshavet
plot(1,1,'-','color',[88/255 122/255 188/255],'linewidth',1.7)
plot(1,1,'-','color',[204/255 102/255 0/255],'linewidth',1.7)
plot(1,1,'-','color',[157/255 14/255 0/255],'linewidth',1.7)
l=legend('1979-1999','1999-2009','2009-2017','location','southeast','AutoUpdate','off'); set(l,'box','on')
h = m_pcolor(lon,lat,Q_mean_refPeriod); set(h, 'EdgeColor', 'none')
cmap = flipud(cbrewer('div','RdBu',30,'linear')); colormap(cmap(11:end,:)), caxis([-1000 3000])
m_grid('xtick',0:20:80,'xticklabels',[],'ytick',70:2:84,'yticklabels',[]); % m_grid('xtick',[],'ytick',[])
m_usercoast('PerfectForHighVertResAtlas_intermediateCoastline','patch',[.6 .6 .6],'edgecolor','none'), m_etopo2('contour',[-200 -1000],'color',[.7 .7 .7])
cb = colorbar('location','southoutside'); xlabel(cb,'Heat content (MJ m^{-2})'), clear cb
% adding contours for june sea ice extent
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19791999,[0.15 0.15],'color',[1 1 1],'linewidth',2.5)
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19791999,[0.15 0.15],'color',[88/255 122/255 188/255],'linewidth',1.7)
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19992009,[0.15 0.15],'color',[1 1 1],'linewidth',2.5)
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19992009,[0.15 0.15],'color',[204/255 102/255 0/255],'linewidth',1.7)
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune20092017,[0.15 0.15],'color',[1 1 1],'linewidth',2.5)
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune20092017,[0.15 0.15],'color',[157/255 14/255 0/255],'linewidth',1.7)
tmpind=indgeo+1-1; tmpind(40,54)=1; tmpind(35,53:54)=1; tmpind(48,48)=0; m_contour(lon,lat,tmpind,[1 1],'color',[.2 .2 .2],'linewidth',1), clear tmpind % legger til kontur for området Study area
m_text(11.5,79.1,'Svalbard'), m_text(41,81.9,{'Franz Josef';'  Land'}), m_text(57,75,{'Novaya';'Zemlya'})

% plot anomalies
% Q anom 2006
figure, hold on, m_proj('Stereographic','lon',38,'lat',76,'rad',[18 68],'rec','on'); % Hele Barentshavet
plot(1,1,'-','color',[88/255 122/255 188/255],'linewidth',1.7)
plot(1,1,'-','color',[204/255 102/255 0/255],'linewidth',1.7)
plot(1,1,'-','color',[157/255 14/255 0/255],'linewidth',1.7)
l=legend('1979-1999','1999-2009','2006','location','southeast','AutoUpdate','off'); set(l,'box','on')
cmap = flipud(cbrewer('div','RdBu',30,'linear')); colormap(cmap(11:end,:)), caxis([-500 1500])
m_grid('xtick',0:20:80,'xticklabels',[],'ytick',70:2:84,'yticklabels',[]); % m_grid('xtick',[],'ytick',[])
m_usercoast('PerfectForHighVertResAtlas_intermediateCoastline','patch',[.6 .6 .6],'edgecolor','none'), 
cb = colorbar('location','southoutside'); xlabel(cb,'Heat content anomaly (MJ m^{-2})'), clear cb
h = m_pcolor(lon,lat,Q_ano2006); set(h, 'EdgeColor', 'none')
m_etopo2('contour',[-200 -1000],'color',[.7 .7 .7])
% adding June SIE contours for the early period, 1979-1999
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19791999,[0.15 0.15],'color',[1 1 1],'linewidth',2.5);
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19791999,[0.15 0.15],'color',[88/255 122/255 188/255],'linewidth',1.7)
% adding June SIE contours for the reference period, 1999-2009
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19992009,[0.15 0.15],'color',[1 1 1],'linewidth',2.5);
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19992009,[0.15 0.15],'color',[204/255 102/255 0/255],'linewidth',1.7)
% adding June SIE contours for the anomalous year, 2006
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune2006,[0.15 0.15],'color',[1 1 1],'linewidth',2.5);
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune2006,[0.15 0.15],'color',[157/255 14/255 0/255],'linewidth',1.7)
tmpind=indgeo+1-1; tmpind(40,54)=1; tmpind(35,53:54)=1; tmpind(48,48)=0; m_contour(lon,lat,tmpind,[1 1],'color',[.2 .2 .2],'linewidth',1), clear tmpind % legger til kontur for området Study area
m_text(11.5,79.1,'Svalbard'), m_text(41,81.9,{'Franz Josef';'  Land'}), m_text(57,75,{'Novaya';'Zemlya'})

% Q anom 2012
figure, hold on, m_proj('Stereographic','lon',38,'lat',76,'rad',[18 68],'rec','on'); % Hele Barentshavet
plot(1,1,'-','color',[88/255 122/255 188/255],'linewidth',1.7)
plot(1,1,'-','color',[204/255 102/255 0/255],'linewidth',1.7)
plot(1,1,'-','color',[157/255 14/255 0/255],'linewidth',1.7)
l=legend('1979-1999','1999-2009','2012','location','southeast','AutoUpdate','off'); set(l,'box','on')
cmap = flipud(cbrewer('div','RdBu',30,'linear')); colormap(cmap(11:end,:)), caxis([-500 1500])
m_grid('xtick',0:20:80,'xticklabels',[],'ytick',70:2:84,'yticklabels',[]); % m_grid('xtick',[],'ytick',[])
m_usercoast('PerfectForHighVertResAtlas_intermediateCoastline','patch',[.6 .6 .6],'edgecolor','none'), 
cb = colorbar('location','southoutside'); xlabel(cb,'Heat content anomaly (MJ m^{-2})'), clear cb
h = m_pcolor(lon,lat,Q_ano2012); set(h, 'EdgeColor', 'none')
m_etopo2('contour',[-200 -1000],'color',[.7 .7 .7])
% adding June SIE contours for the early period, 1979-1999
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19791999,[0.15 0.15],'color',[1 1 1],'linewidth',2.5);
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19791999,[0.15 0.15],'color',[88/255 122/255 188/255],'linewidth',1.7)
% adding June SIE contours for the reference period, 1999-2009
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19992009,[0.15 0.15],'color',[1 1 1],'linewidth',2.5);
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19992009,[0.15 0.15],'color',[204/255 102/255 0/255],'linewidth',1.7)
% adding June SIE contours for the anomalous year, 2012
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune2012,[0.15 0.15],'color',[1 1 1],'linewidth',2.5);
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune2012,[0.15 0.15],'color',[157/255 14/255 0/255],'linewidth',1.7)
tmpind=indgeo+1-1; tmpind(40,54)=1; tmpind(35,53:54)=1; tmpind(48,48)=0; m_contour(lon,lat,tmpind,[1 1],'color',[.2 .2 .2],'linewidth',1), clear tmpind % legger til kontur for området Study area
m_text(11.5,79.1,'Svalbard'), m_text(41,81.9,{'Franz Josef';'  Land'}), m_text(57,75,{'Novaya';'Zemlya'})

% Q anom 2016
figure, hold on, m_proj('Stereographic','lon',38,'lat',76,'rad',[18 68],'rec','on'); % Hele Barentshavet
plot(1,1,'-','color',[88/255 122/255 188/255],'linewidth',1.7)
plot(1,1,'-','color',[204/255 102/255 0/255],'linewidth',1.7)
plot(1,1,'-','color',[157/255 14/255 0/255],'linewidth',1.7)
l=legend('1979-1999','1999-2009','2016','location','southeast','AutoUpdate','off'); set(l,'box','on')
cmap = flipud(cbrewer('div','RdBu',30,'linear')); colormap(cmap(11:end,:)), caxis([-500 1500])
m_grid('xtick',0:20:80,'xticklabels',[],'ytick',70:2:84,'yticklabels',[]); % m_grid('xtick',[],'ytick',[])
m_usercoast('PerfectForHighVertResAtlas_intermediateCoastline','patch',[.6 .6 .6],'edgecolor','none'), 
h = m_pcolor(lon,lat,Q_ano2016); set(h, 'EdgeColor', 'none')
cb = colorbar('location','southoutside'); xlabel(cb,'Heat content anomaly (MJ m^{-2})'), clear cb
m_etopo2('contour',[-200 -1000],'color',[.7 .7 .7])
% adding June SIE contours for the early period, 1979-1999
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19791999,[0.15 0.15],'color',[1 1 1],'linewidth',2.5);
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19791999,[0.15 0.15],'color',[88/255 122/255 188/255],'linewidth',1.7)
% adding June SIE contours for the reference period, 1999-2009
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19992009,[0.15 0.15],'color',[1 1 1],'linewidth',2.5);
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune19992009,[0.15 0.15],'color',[204/255 102/255 0/255],'linewidth',1.7)
% adding June SIE contours for the anomalous year, 2016
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune2016,[0.15 0.15],'color',[1 1 1],'linewidth',2.5);
m_contour(ICEconc_lon,ICEconc_lat,ICEconc_meanJune2016,[0.15 0.15],'color',[157/255 14/255 0/255],'linewidth',1.7)
tmpind=indgeo+1-1; tmpind(40,54)=1; tmpind(35,53:54)=1; tmpind(48,48)=0; h2=m_contour(lon,lat,tmpind,[1 1],'color',[.2 .2 .2],'linewidth',1); clear tmpind % legger til kontur for området Study area
m_text(11.5,79.1,'Svalbard'), m_text(41,81.9,{'Franz Josef';'  Land'}), m_text(57,75,{'Novaya';'Zemlya'})
