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
            load(['D:\projectnew\ionodata5\WWV_ionosphere_','Mar','_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat']);
           
            
            % convert plasma frequency grid to  electron density in electrons/cm^3


            iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
            iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;
            %
        
            % call raytrace
            %
            fr=10;
            bearing_angle=74:0.05:76.95;
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
                    OX_mode, nhops, tol);%X mode rays
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
                %   num = length(ray_N(rayId).lat);
                %   ground_range = zeros(1, num);
                %   lat = ray_N(rayId).lat;
                %   lon = ray_N(rayId).lon;
                %   ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
                %       origin_long,'wgs84')/1000.0;
                %   ray_N(rayId).ground_range = ground_range;
                % end
                for i=1:length(elevs)% select useful ray data and form structure
                    ray_O_final(i).lat=ray_O(i).lat;
                    ray_O_final(i).lon=ray_O(i).lon;
                    ray_O_final(i).height=ray_O(i).height;
                    ray_O_final(i).group_range=ray_O(i).group_range;
                    ray_O_final(i).ground_range=ray_O(i).ground_range;
                    ray_O_final(i).absorption=ray_O(i).absorption;
                    ray_O_final(i).label=ray_data_O(i).ray_label;
                    ray_O_final(i).hopnumber=ray_data_O(i).nhops_attempted;
                    ray_O_final(i).refractive_index=ray_O(i).refractive_index;
                    ray_O_final(i).doppler_shift=ray_data_O(i).Doppler_shift;
                    ray_O_final(i).phase_path=ray_O(i).phase_path;
                    ray_O_final(i).wavenorm_B_angle=ray_O(i).wavenorm_B_angle;
                    ray_O_final(i).TEC_path=ray_data_O(i).TEC_path;
                    ray_O_final(i).electron_density=ray_O(i).electron_density;
                    
                    ray_X_final(i).lat=ray_X(i).lat;
                    ray_X_final(i).lon=ray_X(i).lon;
                    ray_X_final(i).height=ray_X(i).height;
                    ray_X_final(i).group_range=ray_X(i).group_range;
                    ray_X_final(i).ground_range=ray_X(i).ground_range;
                    ray_X_final(i).absorption=ray_X(i).absorption;
                    ray_X_final(i).lable=ray_data_X(i).ray_label;
                    ray_X_final(i).hopnumber=ray_data_X(i).nhops_attempted;
                    ray_X_final(i).refractive_index=ray_X(i).refractive_index;
                    ray_X_final(i).doppler_shift=ray_data_X(i).Doppler_shift;
                    ray_X_final(i).phase_path=ray_X(i).phase_path;
                    ray_X_final(i).wavenorm_B_angle=ray_X(i).wavenorm_B_angle;
                    ray_X_final(i).TEC_path=ray_data_X(i).TEC_path;
                    ray_X_final(i).electron_density=ray_X(i).electron_density;
                end
                filename=['D:\projectnew\ionodata5\raydata_',num2str(fr),'_',num2str(bearing_angle(bearnumber)),'_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'];
                save(filename,'ray_O_final','ray_X_final');%save the raydata

                toc
            end
        end
    end
end
%%
clear
wgs84 = wgs84Ellipsoid;
%bearing_angle=77:0.05:77.95;
bearing_angle=77:0.05:80;
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
                filename=['D:\projectnew\ionodata5\raydata_',num2str(fr),'_',num2str(bearing_angle(bearnumber)),'_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'];
                load(filename);
                i=1;
                m=1;
                for elevsnum=1:length(elevs)
                    
        %             if ray_O_final(elevsnum).ray_label==1
                        value=ray_O_final(elevsnum).height(end);
                        landlocation=find(ray_O_final(elevsnum).height==0);
                        lat_land=ray_O_final(elevsnum).lat(landlocation);
                        lon_land=ray_O_final(elevsnum).lon(landlocation);
                        [bearing,elev,slantRange] = geodetic2aer(lat_land,lon_land,h,rec_lat,rec_lon,h,wgs84);
                      
                        distance=find(slantRange<3000000);
                        if ~isempty(distance)
                            ray_end=landlocation(distance);
                            ray_O_select(i).lat=ray_O_final(elevsnum).lat(1:ray_end);
                            ray_O_select(i).lon=ray_O_final(elevsnum).lon(1:ray_end);
                            ray_O_select(i).height=ray_O_final(elevsnum).height(1:ray_end);
                            ray_O_select(i).group_range=ray_O_final(elevsnum).group_range(1:ray_end);
                            ray_O_select(i).ground_range=ray_O_final(elevsnum).ground_range(1:ray_end);
                            ray_O_select(i).phase_path=ray_O_final(elevsnum).phase_path(1:ray_end);
                            ray_O_select(i).refractive_index=ray_O_final(elevsnum).refractive_index(1:ray_end);
                            ray_O_select(i).elevs=elevs(elevsnum);
                            ray_O_select(i).bearing=bearing_angle(bearnumber);
                            ray_O_select(i).date=date_matrix(timecount);
                            ray_O_select(i).month=month_matrix(timecount);
                            i=i+1;
                        end
        %             end
                    
        %             if ray_X_final(elevsnum).ray_label==1
                        value_X=ray_X_final(elevsnum).height(end);
                        landlocation_X=find(ray_X_final(elevsnum).height==0);
                        lat_land_X=ray_X_final(elevsnum).lat(landlocation_X);
                        lon_land_X=ray_X_final(elevsnum).lon(landlocation_X);
                        [bearing,elev,slantRange_X] = geodetic2aer(lat_land_X,lon_land_X,h,rec_lat,rec_lon,h,wgs84);
                        distance_X=find(slantRange_X<3000000);
                        if ~isempty(distance_X)
                            ray_end_X=landlocation_X(distance_X);
                            ray_X_select(m).lat=ray_X_final(elevsnum).lat(1:ray_end_X);
                            ray_X_select(m).lon=ray_X_final(elevsnum).lon(1:ray_end_X);
                            ray_X_select(m).height=ray_X_final(elevsnum).height(1:ray_end_X);
                            ray_X_select(m).ground_range=ray_X_final(elevsnum).ground_range(1:ray_end_X);
                            ray_X_select(m).group_range=ray_X_final(elevsnum).group_range(1:ray_end_X);
                            ray_X_select(m).phase_path=ray_X_final(elevsnum).phase_path(1:ray_end_X);
                            ray_X_select(m).refractive_index=ray_X_final(elevsnum).refractive_index(1:ray_end_X);
                            ray_X_select(m).elevs=elevs(elevsnum);
                            ray_X_select(m).bearing=bearing_angle(bearnumber);
                            ray_X_select(m).date=date_matrix(timecount);
                            ray_X_select(m).month=month_matrix(timecount);
                            m=m+1;
                        end
                end
            end
            save(['D:\projectnew\ionodata5\received_ray_information','_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'],'ray_O_select','ray_X_select');
        end
%         save(['D:\projectnew\ionodata5\received_ray_information','_',num2str(month),'_',num2str(date),'_',num2str(hour),'_',num2str(minute),'.mat'],'ray_O_select','ray_X_select');
     end
end