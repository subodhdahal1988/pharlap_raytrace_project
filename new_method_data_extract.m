clear
close all
    %parpool(2)
%    parpool(5)   
month_matrix=[3];
date_matrix=[20];
R12_matrix=[52];
minutes_matrix=0:10:59;
hour_matrix=0:1:23;
tic
for timecount=1:length(date_matrix)
    for min=1:length(minutes_matrix)
        for hor=1:length(hour_matrix)
            date=date_matrix(timecount);
            month=month_matrix(timecount);
            minute=minutes_matrix(min);
            hour=hour_matrix(hor);
            load(['E:\projectnew\ionodata5\WWV_ionosphere_','Mar','_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat']);
           
            
            % convert plasma frequency grid to  electron density in electrons/cm^3


            iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
            iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;
            %
        
            % call raytrace
            %
            fr=10;
            bearing_angle=78:0.01:80;
            %bearing_angle=77:0.05:77.95;
            %bearing_angle=78:0.05:80;
            elevs = 3:0.1:90;
            origin_lat = 40+40/60+47.7/3600;             % latitude of the start point of rays
            origin_long = -(105+2/60+27.4/3600);            % longitude of the start point of rays
            origin_ht = 1.583;                % altitude of the start point of rays
            doppler_flag = 0;               % interested in Doppler shift

            %
            % generate ionospheric, geomagnetic and irregularity grids
            %
            ht_start = 90;          % start height for ionospheric grid (km)
            ht_inc = 5;             % height increment step length (km)
            num_ht = 80;
            lat_start = -50;
            lat_inc = 2;
            num_lat = 50;
            lon_start=0;
            lon_inc = 4;
            num_lon = 90;% change here
            iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
                ht_start, ht_inc, num_ht, ];

            B_ht_start = ht_start;          % start height for geomagnetic grid (km)
            B_ht_inc = 5;                  % height increment (km)
            B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
            B_lat_start = lat_start;
            B_lat_inc = 2.0;
            B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
            B_lon_start = lon_start;
            B_lon_inc = 4.0;
            B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);
            geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
                B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

            
            freqs = ones(size(elevs))*fr;   % frequency (MHz)
            tic
            for bearnumber=1:length(bearing_angle)
                
                ray_O_final=struct('lat',[],'lon',[],'height',[],'group_range',[],...
                    'ground_range',[],'absorption',[],'label',[],'hopnumber',[],'refractive_index',[]...
                    ,'phase_path',[],'wavenorm_B_angle',[],'TEC',[],'electron_density',[]);
                repmat(ray_O_final,[1,length(elevs)]);
                ray_X_final=struct('lat',[],'lon',[],'height',[],'group_range',[],...
                    'ground_range',[],'absorption',[],'label',[],'hopnumber',[],...
                    'refractive_index',[],'phase_path',[],'wavenorm_B_angle',[],...
                    'TEC_path',[],'electron_density',[]);
                repmat(ray_X_final,[1,length(elevs)]);
                ray_bears =ones(size(elevs))*bearing_angle(bearnumber);
                nhops = 2;                  % number of hops
                tol = [1e-5 0.01 5];       % ODE solver tolerance and min max stepsizes
                num_elevs = length(elevs);
        
                OX_mode = 1;
                [ray_data_O, ray_O, ray_state_vec_O] = ...
                    raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
                    OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
                    collision_freq, iono_grid_parms, Bx, By, Bz, ...
                    geomag_grid_parms);%O mode rays
                for rayId=1:num_elevs
                    num = length(ray_O(rayId).lat);
                    ground_range = zeros(1, num);
                    lat = ray_O(rayId).lat;
                    lon = ray_O(rayId).lon;
                    ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
                        origin_long,'wgs84')/1000.0;
                    ray_O(rayId).ground_range = ground_range;
                end
        
                OX_mode = -1;
                [ray_data_X, ray_X, ray_sv_X] = ...
                    raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
                    OX_mode, nhops, tol,iono_en_grid, iono_en_grid_5, ...
                    collision_freq, iono_grid_parms, Bx, By, Bz, ...
                    geomag_grid_parms);%X mode rays
                for rayId=1:num_elevs
                    num = length(ray_X(rayId).lat);
                    ground_range = zeros(1, num);
                    lat = ray_X(rayId).lat;
                    lon = ray_X(rayId).lon;
                    ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
                        origin_long,'wgs84')/1000.0;
                    ray_X(rayId).ground_range = ground_range;
                end
                
                
                %      OX_mode = 0;
                % [ray_data_N, ray_N, ray_sv_N] = ...
                %   raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
                %               OX_mode, nhops, tol);
                
                % fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
                %         [ray_data_N.NRT_elapsed_time], NRT_total_time)
                
                % for rayId=1:num_elevs
                  %num = length(ray_N(rayId).lat);
                %   ground_range = zeros(1, num);
                %   lat = ray_N(rayId).lat;
                %   lon = ray_N(rayId).lon;
                %   ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
                %       origin_long,'wgs84')/1000.0;
                %   ray_N(rayId).ground_range = ground_range;
                % end
                for i=1:length(elevs)% select useful ray data and form structure
                    ray_O_final(i).lat=ray_data_O(i).lat;
                    ray_O_final(i).lon=ray_data_O(i).lon;
%                     ray_O_final(i).height=ray_O(i).height;
                    ray_O_final(i).group_range=ray_data_O(i).group_range;
                    ray_O_final(i).ground_range=ray_data_O(i).ground_range;
%                     ray_O_final(i).absorption=ray_O(i).absorption;
                     ray_O_final(i).label=ray_data_O(i).ray_label;
                    ray_O_final(i).hopnumber=ray_data_O(i).nhops_attempted;
%                     ray_O_final(i).refractive_index=ray_O(i).refractive_index;
                    ray_O_final(i).doppler_shift=ray_data_O(i).Doppler_shift;
                    ray_O_final(i).phase_path=ray_data_O(i).phase_path;
%                     ray_O_final(i).wavenorm_B_angle=ray_O(i).wavenorm_B_angle;
                    ray_O_final(i).TEC_path=ray_data_O(i).TEC_path;
%                     ray_O_final(i).electron_density=ray_O(i).electron_density;
                    
                    ray_X_final(i).lat=ray_data_X(i).lat;
                    ray_X_final(i).lon=ray_data_X(i).lon;
%                     ray_X_final(i).height=ray_X(i).height;
                    ray_X_final(i).group_range=ray_data_X(i).group_range;
                    ray_X_final(i).ground_range=ray_data_X(i).ground_range;
%                     ray_X_final(i).absorption=ray_X(i).absorption;
                    ray_X_final(i).label=ray_data_X(i).ray_label;
                    ray_X_final(i).hopnumber=ray_data_X(i).nhops_attempted;
%                     ray_X_final(i).refractive_index=ray_X(i).refractive_index;
                    ray_X_final(i).doppler_shift=ray_data_X(i).Doppler_shift;
                    ray_X_final(i).phase_path=ray_X(i).phase_path;
%                     ray_X_final(i).wavenorm_B_angle=ray_X(i).wavenorm_B_angle;
                    ray_X_final(i).TEC_path=ray_data_X(i).TEC_path;
%                     ray_X_final(i).electron_density=ray_X(i).electron_density;
                end
                filename=['E:\projectnew\pharlap_raytrace_project\local_data\raydata_',num2str(fr),'_',num2str(bearing_angle(bearnumber)),'_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'];
                save(filename,'ray_O_final','ray_X_final');%save the raydata

                toc
            end
        end
    end
end
%% collecting data of ray received 

clear
wgs84 = wgs84Ellipsoid;
%bearing_angle=77:0.05:77.95;
bearing_angle=78:0.01:80;
fr=10;   % frequency (MHz)
rec_lat=40+44/60+30/3600;
rec_lon=-(74+10/60+43/3600);
elevs=3:0.1:90;
h=0;
month_matrix=[3];
date_matrix=[20];
R12_matrix=[52];
minutes_matrix=0:10:59;
hour_matrix=0:1:23;
         ray_O_select=struct('lat',[],'lon',[],'height',[],'group_range',[],...
            'ground_range',[],'absorption',[],'label',[],'hopnumber',[],'refractive_index',[]...
            ,'phase_path',[],'wavenorm_B_angle',[],'TEC',[],'electron_density',[]);
        repmat(ray_O_select,[1,length(elevs)]);
        ray_X_select=struct('lat',[],'lon',[],'height',[],'group_range',[],...
            'ground_range',[],'absorption',[],'label',[],'hopnumber',[],...
            'refractive_index',[],'phase_path',[],'wavenorm_B_angle',[],...
            'TEC_path',[],'electron_density',[]);
        repmat(ray_X_select,[1,length(elevs)]);
        tic 
for timecount=1:length(month_matrix)
     for min=1:length(minutes_matrix)
        for hor=1:length(hour_matrix)
            date=date_matrix(timecount);
            month=month_matrix(timecount);
            minute=minutes_matrix(min);
            hour=hour_matrix(hor);
            for bearnumber=1:length(bearing_angle)
                filename=['E:\projectnew\pharlap_raytrace_project\local_data\raydata_',num2str(fr),'_',num2str(bearing_angle(bearnumber)),'_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'];
                load(filename);
                i=1;
                m=1;
                for elevsnum=1:length(elevs)
                    
                      if ray_O_final(elevsnum).label(end)==1
                            
                            ray_O_select(i).lat=ray_O_final(elevsnum).lat(end);
                            ray_O_select(i).lon=ray_O_final(elevsnum).lon(end);
                            ray_O_select(i).phase_path=ray_O_final(elevsnum).phase_path(end);
    
                            i=i+1;
                      end
                      if ray_X_final(elevsnum).label(end)==1
                          ray_X_select(m).lat=ray_X_final(elevsnum).lat(end);
                          ray_X_select(m).lon=ray_X_final(elevsnum).lon(end);
                          ray_X_select(m).phase_path=ray_X_final(elevsnum).phase_path(end);
                          m=m+1;
                      end
                end
                save(['E:\projectnew\pharlap_raytrace_project\local_data\received_ray_information','_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'],'ray_O_select','ray_X_select');
            end
        end
       
     end
%       save(['D:\projectnew\pharlap_raytrace_project\local_data\received_ray_information','_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'],'ray_O_select','ray_X_select');
end

%                         [bearing,elev,slantRange] = geodetic2aer(lat_land,lon_land,h,rec_lat,rec_lon,h,wgs84);
%                         
%                       
%                         distance=find(slantRange<30000);
%                         if ~isempty(distance)
%                             ray_end=landlocation(distance);
%                             ray_O_select(i).lat=ray_O_final(elevsnum).lat(1:ray_end);
%                             ray_O_select(i).lon=ray_O_final(elevsnum).lon(1:ray_end);
% %                             ray_O_select(i).height=ray_O_final(elevsnum).height(1:ray_end);
%                             ray_O_select(i).group_range=ray_O_final(elevsnum).group_range(1:ray_end);
%                             ray_O_select(i).ground_range=ray_O_final(elevsnum).ground_range(1:ray_end);
%                             ray_O_select(i).phase_path=ray_O_final(elevsnum).phase_path(1:ray_end);
% %                             ray_O_select(i).refractive_index=ray_O_final(elevsnum).refractive_index(1:ray_end);
%                             ray_O_select(i).elevs=elevs(elevsnum);
%                             ray_O_select(i).bearing=bearing_angle(bearnumber);
%                             ray_O_select(i).date=date_matrix(timecount);
%                             ray_O_select(i).month=month_matrix(timecount);
%                             i=i+1;
%                         end
%         %             end
%                     
%                     if ray_X_final(elevsnum).label==1
%                         
%                         lat_land_X=ray_X_final(elevsnum).lat(end);
%                         lon_land_X=ray_X_final(elevsnum).lon(end);
%                         [bearing,elev,slantRange_X] = geodetic2aer(lat_land_X,lon_land_X,h,rec_lat,rec_lon,h,wgs84);
%                         distance_X=find(slantRange_X<30000);
%                         if ~isempty(distance_X)
%                             ray_end_X=landlocation_X(distance_X);
%                             ray_X_select(m).lat=ray_X_final(elevsnum).lat(1:ray_end_X);
%                             ray_X_select(m).lon=ray_X_final(elevsnum).lon(1:ray_end_X);
% %                             ray_X_select(m).height=ray_X_final(elevsnum).height(1:ray_end_X);
%                             ray_X_select(m).ground_range=ray_X_final(elevsnum).ground_range(1:ray_end_X);
%                             ray_X_select(m).group_range=ray_X_final(elevsnum).group_range(1:ray_end_X);
%                             ray_X_select(m).phase_path=ray_X_final(elevsnum).phase_path(1:ray_end_X);
% %                             ray_X_select(m).refractive_index=ray_X_final(elevsnum).refractive_index(1:ray_end_X);
%                             ray_X_select(m).elevs=elevs(elevsnum);
%                             ray_X_select(m).bearing=bearing_angle(bearnumber);
%                             ray_X_select(m).date=date_matrix(timecount);
%                             ray_X_select(m).month=month_matrix(timecount);
%                             m=m+1;
% 
%                         end
%                     end
%                   end
%                 end
%             end
%         end
%      end
%      save(['D:\projectnew\pharlap_raytrace_project\local_data\received_ray_information','_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'],'ray_O_select','ray_X_select');
% end
% 
%                         
%                
%             
%            
%      
%    
%% finding slant range of ray landed from NJIT
clear
lat0=40.7417;
lon0=-74.1786;
%lat0=40+44/60+30/3600;
%lon0=-(74+10/60+43/3600);
h=0;

wgs84 = wgs84Ellipsoid;
month_matrix=[3];
date_matrix=[20];
minute_mat=0:10:59;
hor_mat=0:1:23;

ray_O_received_lowestdistance=struct('lat',[],'lon',[],'ground_range',[],'phase_path',[],'slantrange',[]);
ray_X_received_lowestdistance=struct('lat',[],'lon',[],'ground_range',[],'phase_path',[],'slantrange',[]);
i=1;
m=1;
for timecount=1:length(month_matrix)
    for l=1:length(hor_mat)
        for j=1:length(minute_mat)
            date=date_matrix(timecount);
            month=month_matrix(timecount);
            hour=hor_mat(l);
            minute=minute_mat(j);
            filename=['E:\projectnew\pharlap_raytrace_project\local_data\received_ray_information','_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'];
            load(filename);
            
            for num1=1:length(ray_O_select)
                
                    
                    lat_land=ray_O_select(num1).lat;
                    lon_land=ray_O_select(num1).lon;
                    
                    [az,elev,slantRange_O]=geodetic2aer(lat_land,lon_land,h,lat0,lon0,h,wgs84);

                    ray_O_received_lowestdistance(num1).lat=ray_O_select(num1).lat;
                    ray_O_received_lowestdistance(num1).lon=ray_O_select(num1).lon;

                    ray_O_received_lowestdistance(num1).phase_path=ray_O_select(num1).phase_path;
                    ray_O_received_lowestdistance(num1).slantrange=slantRange_O;

%                     ray_O_received_lowestdistance(i).month=month;
%                     ray_O_received_lowestdistance(i).date=date;
%                     ray_O_received_lowestdistance(i).hour=hour;
%                     ray_O_received_lowestdistance(i).minute=minute;
%                     
                    i=i+1;
        
            end

   
           
            

             for num2=1:length(ray_X_select)
                   
                    lat_land_X=ray_X_select(num2).lat;
                    lon_land_X=ray_X_select(num2).lon;
                    [az,elev,slantRange_X]=geodetic2aer(lat_land_X,lon_land_X,h,lat0,lon0,h,wgs84);
                    


                    

                    ray_X_received_lowestdistance(num2).lat=ray_X_select(num2).lat;
                    ray_X_received_lowestdistance(num2).lon=ray_X_select(num2).lon;
                    
                    ray_X_received_lowestdistance(num2).phase_path=ray_X_select(num2).phase_path;
                    ray_X_received_lowestdistance(num2).slantrange=slantRange_X;
                    
%                     ray_X_received_lowestdistance(m).month=month;
%                     ray_X_received_lowestdistance(m).date=date;
%                     ray_X_received_lowestdistance(m).hour=hour;
%                     ray_X_received_lowestdistance(m).minute=minute;
                    m=m+1;
             end
             save(['E:\projectnew\pharlap_raytrace_project\local_data\received_ray_closetoreciever','_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'],'ray_X_received_lowestdistance','ray_O_received_lowestdistance')

                    
        end
% %         save(['D:\projectnew\pharlap_raytrace_project\local_data\received_ray_closetoreciever','_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'],'ray_X_received_lowestdistance','ray_O_received_lowestdistance')
       
    end
    
end
%% finding O and X ray landed near NJIT
clear
lat0=40.7417;
lon0=-74.1786;
%lat0=40+44/60+30/3600;
%lon0=-(74+10/60+43/3600);
h=0;

wgs84 = wgs84Ellipsoid;
month_matrix=[3];
date_matrix=[20];
minute_mat=0:10:59;
hor_mat=0:1:23;

ray_O_received_njit=struct('lat',[],'lon',[],'ground_range',[],'phase_path',[],'slantrange',[]);
ray_X_received_njit=struct('lat',[],'lon',[],'ground_range',[],'phase_path',[],'slantrange',[]);
i=1;
m=1;
for timecount=1:length(month_matrix)
    for l=1:length(hor_mat)
        for j=1:length(minute_mat)
            date=date_matrix(timecount);
            month=month_matrix(timecount);
            hour=hor_mat(l);
            minute=minute_mat(j);
            filename=['E:\projectnew\pharlap_raytrace_project\local_data\received_ray_closetoreciever','_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'];
            load(filename);


            for num1=1:length(ray_O_received_lowestdistance)
                    n1(num1)=ray_O_received_lowestdistance(num1).slantrange;
            end
                    position=find(n1==min(n1));
            
                    ray_O_received_njit(i).lat=ray_O_received_lowestdistance(position).lat;
                    ray_O_received_njit(i).lon=ray_O_received_lowestdistance(position).lon;

                    ray_O_received_njit(i).phase_path=ray_O_received_lowestdistance(position).phase_path;
                   

                    ray_O_received_njit(i).month=month;
                    ray_O_received_njit(i).date=date;
                    ray_O_received_njit(i).hour=hour;
                    ray_O_received_njit(i).minute=minute;
%                     
                    i=i+1;
        
            


             for num2=1:length(ray_X_received_lowestdistance)
                   
                    n2(num2)=ray_X_received_lowestdistance(num2).slantrange;
             end
                    
                    position_x=find(n2==min(n2));
                  

                    ray_X_received_njit(m).lat=ray_X_received_lowestdistance(position_x).lat;
                    ray_X_received_njit(m).lon=ray_X_received_lowestdistance(position_x).lon;
                    
                    ray_X_received_njit(m).phase_path=ray_X_received_lowestdistance(position_x).phase_path;
                   
                    
                    ray_X_received_njit(m).month=month;
                    ray_X_received_njit(m).date=date;
                    ray_X_received_njit(m).hour=hour;
                    ray_X_received_njit(m).minute=minute;
                    m=m+1;
        end
    end
    save(['E:\projectnew\pharlap_raytrace_project\local_data\landedray_njit','_',num2str(month),'_',num2str(date),'.mat'],'ray_X_received_njit','ray_O_received_njit')
end
 %% time calculation
 clear
 load('E:\projectnew\pharlap_raytrace_project\local_data\landedray_njit_3_20.mat')
s=1;
for i1=1:length(ray_O_received_njit)
        time_year(i1)=2018;
        time_month(i1)=ray_O_received_njit(i1).month(1);
        time_day(i1)=ray_O_received_njit(i1).date(1);
        time_hour(i1)=ray_O_received_njit(i1).hour(1);
        time_minute(i1)=ray_O_received_njit(i1).minute(1);
        time_second(i1)=00;
   
        now=[time_year(i1),time_month(i1),time_day(i1),time_hour(i1),time_minute(i1),time_second(i1)];
        ray_O_received_njit(i1).time=now;
        %UT1(i)=datetime(ray_O_received_lowestdistance(i).time);
         ray_O_received_njit(s).UT=datetime(ray_O_received_njit(i1).time);
         ray_X_received_njit(s).UT=datetime(ray_O_received_njit(i1).time);
         s=s+1;
end
                   

 %%
%lat vs lon plot (for O ray)

for e=1:length(ray_X_received_njit)
    
    timo(e)=ray_O_received_njit(e).UT;
    O_lat(e)=ray_O_received_njit(e).lat;
    O_lot(e)=ray_O_received_njit(e).lon;

    yyaxis left
    plot(timo,O_lat,'color','b','Marker','o');
    ylabel('latitude')

    hold on
    yyaxis right                
    plot(timo,O_lot,'color','r','Marker','x');




    e=e+1;
    xlabel('UT')
    ylabel('longitude')
    

    title('latitude vs longitude for O-rays landed closest NJIT(2018-03-20)')
    
end
   

             
 %%
%lat vs lon plot (for X ray)

for e=1:length(ray_X_received_njit)
    
    timx(e)=ray_X_received_njit(e).UT;
    X_lat(e)=ray_X_received_njit(e).lat;
    X_lot(e)=ray_X_received_njit(e).lon;

    yyaxis left
    plot(timx,X_lat,'color','b','Marker','o');
    ylabel('latitude')

    hold on
    yyaxis right                
    plot(timx,X_lot,'color','r','Marker','x');




    e=e+1;
    xlabel('UT')
    ylabel('longitude')
    

    title('latitude vs longitude for X-rays landed closest NJIT(2018-03-20)')
    
end
%% geographic plot
lat0=40.7417;
lon0=-74.1786;
geoplot(X_lat,X_lot,"om",MarkerFaceColor="auto",Marker="*",MarkerSize=15)
hold on
geoplot(O_lat,O_lot,"om",MarkerFaceColor="auto",Marker=".",MarkerSize=15)
hold on
geoplot(lat0,lon0,"om",Marker="diamond",MarkerSize=20,MarkerFaceColor="auto")

legend("x-ray","o-ray","NJIT")
%%
%phase path plot


for pc=1:length(ray_X_received_njit);
    
    tim(pc)=ray_O_received_njit(pc).UT;
    phasepath_0(pc)=ray_O_received_njit(pc).phase_path;
    phasepath_1(pc)=ray_X_received_njit(pc).phase_path;


    
end
%%
%filter plot

    order = 3;
    framelen = 21;

    filt_O=sgolayfilt(phasepath_0,order,framelen);
    filt_X=sgolayfilt(phasepath_1,order,framelen);


    plot(tim,phasepath_0,'color','r','Marker','.');

    hold on
    plot(tim,filt_O,'color','black','Marker','+','LineWidth',2);
    hold on

    plot(tim,phasepath_1,'color','b','Marker','x');

    hold on
    plot(tim,filt_X,'color','cyan','Marker','^','LineWidth',2);



  
    


      xlabel('UT')
      ylabel('Phase path(Km)')



      ax = gca;
      ax.XGrid = 'on';
      ax.YGrid = 'off';

    %       %xlim(0,25)
    %  
      title('phase path of X and O Rays received at NJIT(2018-03-20)')
     legend('O-ray','fit_O','X-ray','fit_X','Location','bestoutside')
%%
        
    
   
