%% getting closest ray recieved at NJIT 
clear
lat0=40.7423;
lon0=-74.1793;
%lat0=40+44/60+30/3600;
%lon0=-(74+10/60+43/3600);
h=0;

wgs84 = wgs84Ellipsoid;
month_matrix=[3];
date_matrix=[20];
minute_mat=0:10:59;
hor_mat=0:1:23;

ray_O_received_lowestdistance=struct('lat',[],'lon',[],'ground_range',[],'phase_path',[],'height',[]);
ray_X_received_lowestdistance=struct('lat',[],'lon',[],'ground_range',[],'phase_path',[],'height',[]);
i=1;
m=1;
for timecount=1:length(month_matrix)
    for l=1:length(hor_mat)
        for j=1:length(minute_mat)
            date=date_matrix(timecount);
            month=month_matrix(timecount);
            hor=hor_mat(l);
            minute=minute_mat(j);
            filename=['D:\projectnew\ionodata5\received_ray_information','_',num2str(month),'_',num2str(date),'_',num2str(hor),'_',num2str(minute),'.mat'];
            load(filename);
            
            for num1=1:length(ray_O_select)
                
                    landlocation=find(ray_O_select(num1).height==0);
                    lat_land=ray_O_select(num1).lat(landlocation(end));
                    lon_land=ray_O_select(num1).lon(landlocation(end));
                    
                    [az,elev,slantRange]=geodetic2aer(lat_land,lon_land,h,lat0,lon0,h,wgs84);
                    
                    D(num1)=slantRange;

                    
                   
            end
             distance=find(min(D));       
            
                      
                    %ray_end=landlocation(distance);
                    ray_O_received_lowestdistance(i).lat=ray_O_select(distance).lat;
                    ray_O_received_lowestdistance(i).lon=ray_O_select(distance).lon;
                    ray_O_received_lowestdistance(i).ground_range=ray_O_select(distance).ground_range;
                    ray_O_received_lowestdistance(i).phase_path=ray_O_select(distance).phase_path;
                    ray_O_received_lowestdistance(i).height=ray_O_select(distance).height;
                    ray_O_received_lowestdistance(i).month=month;
                    ray_O_received_lowestdistance(i).date=date;
                    ray_O_received_lowestdistance(i).hour=hor;
                    ray_O_received_lowestdistance(i).minute=minute;
                    
                    i=i+1;

   
           
%                

               for num2=1:length(ray_X_select)
                   landlocation_X=find(ray_X_select(num2).height==0);
                    lat_land_X=ray_X_select(num2).lat(landlocation_X(end));
                    lon_land_X=ray_X_select(num2).lon(landlocation_X(end));
                    [slantRange_X]=geodetic2aer(lat_land_X,lon_land_X,h,lat0,lon0,h,wgs84);
                     N(num2)=slantRange_X;
               end
                    distance_X=find(min(N));


                    

                    ray_X_received_lowestdistance(m).lat=ray_X_select(distance_X).lat;
                    ray_X_received_lowestdistance(m).lon=ray_X_select(distance_X).lon;
                    ray_X_received_lowestdistance(m).ground_range=ray_X_select(distance_X).ground_range;
                    ray_X_received_lowestdistance(m).phase_path=ray_X_select(distance_X).phase_path;
                    ray_X_received_lowestdistance(m).height=ray_X_select(distance_X).height;
                    ray_X_received_lowestdistance(m).month=month;
                    ray_X_received_lowestdistance(m).date=date;
                    ray_X_received_lowestdistance(m).hour=hor;
                    ray_X_received_lowestdistance(m).minute=minute;
                    m=m+1;
                    
        end
       
    end
    save(['D:\projectnew\ionodata5\received_ray_closetoreciever','_',num2str(month),'_',num2str(date),'.mat'],'ray_X_received_lowestdistance','ray_O_received_lowestdistance')
end
%save(['D:\projectnew\ionodata5\received_ray_closetoreciever','_',num2str(month),'_',num2str(date),'.mat'],'ray_X_received_lowestdistance','ray_O_received_lowestdistance')
%%
%3D plot  
clear

load('/Volumes/subodh/iondata/received_ray_closetoreciever_3_20.mat')
for k=1:144
        x1=ray_O_received_lowestdistance(k).lat;
        x2=ray_X_received_lowestdistance(k).lat;
        y1=ray_O_received_lowestdistance(k).lon;
        y2=ray_X_received_lowestdistance(k).lon;

        %[xx,yy]=meshgrid(x,y);
        O=ray_O_received_lowestdistance(k).phase_path();
        X=ray_X_received_lowestdistance(k).phase_path();
    
        %surf(xx,yy,z);
        plot3(x1,y1,O,"Color",'r');
        hold on
        plot3(x2,y2,X,"color", 'b');
        xlabel('Latitude')
        ylabel('Longitude')
        zlabel('phase path')
        title('X and O Rays received at NJIT(2018-03-20)')
        legend('O-ray','X-ray','Location','bestoutside')
end


%%
% ray received closest to NJIT
clear
lat0=40.7423;
lon0=-74.1793;
% lat0=40+44/60+30/3600;
% lon0=-(74+10/60+43/3600);
load('D:\projectnew\ionodata5\received_ray_closetoreciever_3_20.mat')

j=1;
for k=1:length(ray_X_received_lowestdistance)
   len_O=length(ray_O_received_lowestdistance(k).lat());
   ray_O_received_lowestdistance(k).length_O=len_O;
   len_X=length(ray_X_received_lowestdistance(k).lon());
   ray_X_received_lowestdistance(k).length_X=len_X;
   
   j=j+1;
end 
 
i=1;
for m=1:length(ray_X_received_lowestdistance)
    len=ray_O_received_lowestdistance(m).length_O;
    len1=ray_X_received_lowestdistance(m).length_X;
    for b=1:len
        lat_O=ray_O_received_lowestdistance(m).lat(b);
        lon_O=ray_O_received_lowestdistance(m).lon(b);
%         x=[lat0;lon0];
%         y=[lat_O;lon_O];
%         d1= norm(x-y);
        d1=sqrt((lat0-lat_O)^2+(lon0-lon_O)^2);
        ray_O_received_lowestdistance(m).distance(b)=d1;
    end
    for c=1:len1
        lat_X=ray_X_received_lowestdistance(m).lat(c);
        lon_X=ray_X_received_lowestdistance(m).lon(c);
%         x1=[lat0;lon0];
%         y1=[lat_X;lon_X];
%         d2= norm(x1-y1);
        d2=sqrt((lat0-lat_X)^2+(lon0-lon_X)^2);
        ray_X_received_lowestdistance(m).distance_X(c)=d2;
    end

    i=i+1;
    m=m+1;     
end

%%
v=1;
for m=1:length(ray_X_received_lowestdistance)
 dis=ray_O_received_lowestdistance(m).distance();
 
 
 pos=find(dis == min(dis));
 
 ray_O_received_lowestdistance(m).phasepath=ray_O_received_lowestdistance(m).phase_path(pos);
 ray_O_received_lowestdistance(m).latitude=ray_O_received_lowestdistance(m).lat(pos);
 ray_O_received_lowestdistance(m).longitude=ray_O_received_lowestdistance(m).lon(pos);


 v=v+1;
end
%%
v1=1;
for w=1:length(ray_X_received_lowestdistance)

 dis1=ray_X_received_lowestdistance(w).distance_X();
 
 
 pos1=find(dis1 == min(dis1));


 ray_X_received_lowestdistance(w).phasepath1=ray_X_received_lowestdistance(w).phase_path(pos1);
 ray_X_received_lowestdistance(w).latitude1=ray_X_received_lowestdistance(w).lat(pos1);
 ray_X_received_lowestdistance(w).longitude1=ray_X_received_lowestdistance(w).lon(pos1);
 v1=v1+1;
end


%%
%time calculation
s=1;
for i1=1:length(ray_X_received_lowestdistance)
        time_year(i1)=2018;
        time_month(i1)=ray_O_received_lowestdistance(i1).month(1);
        time_day(i1)=ray_O_received_lowestdistance(i1).date(1);
        time_hour(i1)=ray_O_received_lowestdistance(i1).hour(1);
        time_minute(i1)=ray_O_received_lowestdistance(i1).minute(1);
        time_second(i1)=00;
   
        now=[time_year(i1),time_month(i1),time_day(i1),time_hour(i1),time_minute(i1),time_second(i1)];
        ray_O_received_lowestdistance(i1).time=now;
        %UT1(i)=datetime(ray_O_received_lowestdistance(i).time);
         ray_O_received_lowestdistance(s).UT=datetime(ray_O_received_lowestdistance(i1).time);
         ray_X_received_lowestdistance(s).UT=datetime(ray_O_received_lowestdistance(i1).time);
         s=s+1;
end
                   

%%
%phase path plot
data_new=struct('time',[],'phaseo',[],'phasex',[]);

for pc=1:length(ray_X_received_lowestdistance);
    
    tim(pc)=ray_O_received_lowestdistance(pc).UT;
    phasepath_0(pc)=ray_O_received_lowestdistance(pc).phasepath;
    phasepath_1(pc)=ray_X_received_lowestdistance(pc).phasepath1;


  
    data_new.time=tim;
    data_new.phaseo=phasepath_0;
    data_new.phasex=phasepath_1;

    
end
%%
%filter plot

    order = 3;
    framelen = 5;

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

%phase_O=data_new.phaseo;
datamatrix=zeros(length(ray_X_received_lowestdistance),10);
for q=1:length(ray_X_received_lowestdistance);
    T=data_new.time(q);
    phasex=data_new.phasex(q);
    phaseo=data_new.phaseo(q);
    
  
    datamatrix(q,1)=datenum(T);
    datamatrix(q,2)=phasex;
    datamatrix(q,3)=phaseo;
end

% for w=1:143
%     [year,month,day,hour,minutes,seconds]=datevec(datamatrix(:,1));
%     datamatrix(w,4)=year(w,1);
%     datamatrix(w,5)=month(w,1);
%     datamatrix(w,6)=day(w,1);
%     datamatrix(w,7)=hour(w,1);
%     datamatrix(w,8)=minutes(w,1);
%     datamatrix(w,9)=seconds(w,1);
%     datamatrix(w,10)=(hour(w,1)+((minutes(w,1)/60))+((seconds(w,1)/3600)));
%     datamatrix(144,10)=24;
% 
% end

%plot(T,phase_O)
figure
num1=fit(datamatrix(:,1),datamatrix(:,2),'poly3','Normalize','on','Robust','Bisquare');

num2=fit(datamatrix(:,1),datamatrix(:,3),'poly3','Normalize','on','Robust','Bisquare');


plot(datamatrix(:,1),datamatrix(:,2),'color','r','Marker','o');
hold on
plot(datamatrix(:,1),datamatrix(:,3),'color','b','Marker','x');
hold on
Oplot=plot(num1,datamatrix(:,1),datamatrix(:,2));
set(Oplot,'color','r');
hold on
xplot=plot(num2,datamatrix(:,1),datamatrix(:,3));
set(xplot,'color','b');
hold off


legend('O-ray','X-ray')


%% phase path curve with fit data with time axis UT
for nm=1:length(ray_X_received_lowestdistance)
    fitdata=num1(datamatrix(1,1):0.0069:datamatrix(144,1));
    fitdata1=num2(datamatrix(1,1):0.0069:datamatrix(144,1));
    
    ray_O_received_lowestdistance(nm).fitdata_O=fitdata(nm);
    ray_X_received_lowestdistance(nm).fitdata_x=fitdata1(nm);

    plot(tim,phasepath_0,'color','r','Marker','o');
    
    hold on
    
    plot(tim,fitdata,'color','r')

    hold on
    
    plot(tim,phasepath_1,'color','b','Marker','x');

    hold on
   
    plot(tim,fitdata1,'color','b')

    legend('O-ray','O-bestfit','X-ray','X-bestfit')

     xlabel('UT')
      ylabel('Phase path(Km)')
      title('Phasepath of rays landed closest to NJIT(2018-03-20)')
end


%%
%lat vs lon plot (for O ray)
for e=1:length(ray_X_received_lowestdistance)
    
    timo(e)=ray_O_received_lowestdistance(e).UT;
    O_lat(e)=ray_O_received_lowestdistance(e).latitude;
    O_lot(e)=ray_O_received_lowestdistance(e).longitude;

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
%geographic position of O-ray
geoplot(X_lat,X_lot,"om",MarkerFaceColor="auto",Marker="*",MarkerSize=15)
hold on
geoplot(O_lat,O_lot,"om",MarkerFaceColor="auto",Marker=".",MarkerSize=15)
hold on
geoplot(lat0,lon0,"om",Marker="diamond",MarkerSize=20,MarkerFaceColor="auto")

legend("x-ray","o-ray","NJIT")


%%
%lat vs lot (x ray)
for bg=1:length(ray_X_received_lowestdistance)
    
    timx(bg)=ray_X_received_lowestdistance(bg).UT;
    X_lat(bg)=ray_X_received_lowestdistance(bg).latitude1;
    X_lot(bg)=ray_X_received_lowestdistance(bg).longitude1;

    yyaxis left
    plot(timx,X_lat,'color','b','Marker','o');
    ylabel('latitude')

    hold on
    yyaxis right                
    plot(timx,X_lot,'color','r','Marker','x');

    
    xlabel('UT')
    ylabel('longitude')
    %xlim(0,24)

    title('latitude vs longitude for X-rays landed closest to NJIT(2018-03-20)')
    
end

