function target = BMinform(N,selected_target_inform)

target = struct('targetid',0,'BMtype',0,'launchSiteid',0,'launchSite',0,'destinationid',0,'destination',0,...
                'launchDirec',0,'launchinform',0,'Cd',0,'time',0,'Alt',0,'Speed',0,'geodetic',[],'x',0) ;

    % launch angle information-------------------------------------------------
    %--------------------------------------------------------------------------
    site1BM1 = [143.7569         0         0         0         0         0;
                 86.0000         0         0         0         0         0]; 
    site1BM2 = [137.6000  146.0000  130.0000  138.0000  158.0000         0;
                 89.3318   88.9787   89.1361   88.5332   88.5886         0];
    site1BM3 = [136.0000  144.0000  128.0000  138.0000  156.0000  140.0000;
                 89.4231   89.1439   89.2564   88.8846   88.8932   88.6330];
    site1BM4 = [112.0000  128.0000  110.0000  124.0000  146.0000  128.0000;
                 89.8514   89.7853   89.8061   89.7102   89.7418   89.6603];
    %-------------------------------------------------------------------------- 
    site2BM2 = [179.0000  176.0000  168.0000  166.0000         0         0;  
                 89.1002   88.6224   89.0156   87.4776         0         0];
    site2BM3 = [178.5000  174.0000  166.0000  164.0000  180.0000  162.0000;
                 89.2543   88.9285   89.1849   88.7251   88.5021   88.4073];
    site2BM4 = [170.0000  168.0000  154.0000  156.0000  176.0000  164.0000;
                 89.8472   89.7646   89.8202   89.7134   89.6836   87.5357];
    %--------------------------------------------------------------------------
    site3BM2 = [200.5000  192.0000  186.0000  178.0000         0         0;
                 89.1772   88.7872   89.1746   88.6089         0         0];
    site3BM3 = [200.0000  190.0000  186.0000  178.0000  192.0000  174.0000;
                 89.3153   89.0410   89.3125   88.9231   88.6322   88.6879];
    site3BM4 = [202.0000  188.0000  188.0000  172.0000  192.0000  170.0000;
                 89.8682   89.7973   87.2690   89.7647   89.7086   89.7147];
    %--------------------------------------------------------------------------
    site4BM2 = [209.2010  199.0000  197.0000  185.0000         0         0;
                 89.1000   88.7188   89.1319   88.5937         0         0];              
    site4BM3 = [209.0000  197.7070  195.0000  185.0000  199.0000  181.0000;
                 89.2576   89.0000   89.2883   88.9152   88.5542   86.3903];
    site4BM4 = [215.0000  199.0000  195.0000  181.0000  199.0000  175.0000;
                 89.8502   89.7866   89.8606   89.7661   89.6976   89.7173];
    %-------------------------------------------------------------------------- 

    %load Cd(drag coefficient)-------------------------------------------------
    %--------------------------------------------------------------------------
    NodongA = importdata('NodongA.mat')                 ;
    NodongAwarhead = importdata('NodongAwarhead.mat')   ;
    ScudB = importdata('ScudB.mat')                     ;
    ScudBwarhead = importdata('ScudBwarhead.mat')       ;
    ScudC = importdata('ScudC.mat')                     ;
    ScudCwarhead = importdata('ScudCwarhead.mat')       ;
    ScudD = importdata('ScudD.mat')                     ;
    ScudDwarhead = importdata('ScudDwarhead.mat')       ;
    alt = [-1000 ;     0;  5000; 10000; 15000; 20000; 25000; 30000; 35000; 40000;...
            45000; 50000; 55000; 60000; 65000; 70000; 75000; 80000; 85000]                  ;
    machNodong = [0 0.5 1 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5 15.5];   
    machD = [0 0.5 1 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5]                    ;
    machC = [0 0.5 1 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5]                                   ;
    machB = [0 0.5 1 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]                                       ;
    %--------------------------------------------------------------------------

for id = 1:1:N
    
    target(id).targetid      = id         ;    
    target(id).BMtype        = selected_target_inform(1,id);
    target(id).launchSiteid  = selected_target_inform(2,id);
    target(id).destinationid = selected_target_inform(3,id);

    switch target(id).BMtype  % BM type1~4 : Scud-B , Scud-C , Scud-D , Nodong-A
        case 1
            target(id).launchinform= struct('GrWeiStg1',5000,'GrWeiStg2',900,'FuWeiStg1',3600 ,'FuWeiStg2',0,...
                                            'ISPstg1',250,'ISPstg2',0,'BurTimStg1',60,'BurTimStg2',1)      ;
            target(id).Cd = struct('diameter',[0.89 0.885],'Cd1',ScudB,'Cd2',ScudBwarhead,'alt',alt,'mach',machB) ;
       
            switch target(id).launchSiteid
                case 1
                    target(id).launchSite= struct('LatH','N','LatDeg',40,'LonH','E','LonDeg',124.82)     ; % Tongchang Dong
                    TempLaunchAngle_inform = site1BM1 ;
                otherwise
                    disp('there are some problems!!(Launch Site)')        
            end
        case 2
            target(id).launchinform= struct('GrWeiStg1',5350,'GrWeiStg2',750,'FuWeiStg1',4200,'FuWeiStg2',0,...
                                            'ISPstg1',270,'ISPstg2',0,'BurTimStg1',87,'BurTimStg2',1)      ;
            target(id).Cd = struct('diameter',[0.9 0.89],'Cd1',ScudC,'Cd2',ScudCwarhead,'alt',alt,'mach',machC) ;  
               
            switch target(id).launchSiteid
                case 1
                    target(id).launchSite= struct('LatH','N','LatDeg',40,'LonH','E','LonDeg',124.82)     ; % Tongchang Dong
                    TempLaunchAngle_inform = site1BM2 ;
                case 2 
                    target(id).launchSite= struct('LatH','N','LatDeg',41.4,'LonH','E','LonDeg',127.12)   ; % Yong jo Ri Site
                    TempLaunchAngle_inform = site2BM2 ;
                case 3
                    target(id).launchSite= struct('LatH','N','LatDeg',40.95,'LonH','E','LonDeg',128.62)  ; % Sangnam Ri
                    TempLaunchAngle_inform = site3BM2 ;
                case 4
                    target(id).launchSite= struct('LatH','N','LatDeg',40.95,'LonH','E','LonDeg',129.32)  ; % Musudan Ri
                    TempLaunchAngle_inform = site4BM2 ;
                otherwise
                    disp('there are some problems!!(Launch Site)')        
            end
        case 3
            % ¿©±â ¹Ù²Þ BurTimStg1 89 -> 891
            target(id).launchinform= struct('GrWeiStg1',5900,'GrWeiStg2',500,'FuWeiStg1',4600,'FuWeiStg2',0,...
                                            'ISPstg1',270,'ISPstg2',0,'BurTimStg1',89,'BurTimStg2',1)      ;
            target(id).Cd = struct('diameter',[0.89 0.885],'Cd1',ScudD,'Cd2',ScudDwarhead,'alt',alt,'mach',machD) ;

            switch target(id).launchSiteid
                case 1
                    target(id).launchSite= struct('LatH','N','LatDeg',40,'LonH','E','LonDeg',124.82)     ; % Tongchang Dong
                    TempLaunchAngle_inform = site1BM3 ;
                case 2 
                    target(id).launchSite= struct('LatH','N','LatDeg',41.4,'LonH','E','LonDeg',127.12)   ; % Yong jo Ri Site
                    TempLaunchAngle_inform = site2BM3 ;
                case 3
                    target(id).launchSite= struct('LatH','N','LatDeg',40.95,'LonH','E','LonDeg',128.62)  ; % Sangnam Ri
                    TempLaunchAngle_inform = site3BM3 ;
                case 4
                    target(id).launchSite= struct('LatH','N','LatDeg',40.95,'LonH','E','LonDeg',129.32)  ; % Musudan Ri
                    TempLaunchAngle_inform = site4BM3 ;
                otherwise
                    disp('there are some problems!!(Launch Site)')        
            end
        case 4
            target(id).launchinform= struct('GrWeiStg1',15200,'GrWeiStg2',800,'FuWeiStg1',13000,'FuWeiStg2',0,...
                                            'ISPstg1',264,'ISPstg2',0,'BurTimStg1',110,'BurTimStg2',1)      ;
            target(id).Cd = struct('diameter',[1.3 1.25],'Cd1',NodongA,'Cd2',NodongAwarhead,'alt',alt,'mach',machNodong) ;

            switch target(id).launchSiteid
                case 1
                    target(id).launchSite= struct('LatH','N','LatDeg',40,'LonH','E','LonDeg',124.82)     ; % Tongchang Dong
                    TempLaunchAngle_inform = site1BM4 ;
                case 2 
                    target(id).launchSite= struct('LatH','N','LatDeg',41.4,'LonH','E','LonDeg',127.12)   ; % Yong jo Ri Site
                    TempLaunchAngle_inform = site2BM4 ;
                case 3
                    target(id).launchSite= struct('LatH','N','LatDeg',40.95,'LonH','E','LonDeg',128.62)  ; % Sangnam Ri
                    TempLaunchAngle_inform = site3BM4 ;
                case 4
                    target(id).launchSite= struct('LatH','N','LatDeg',40.95,'LonH','E','LonDeg',129.32)  ; % Musudan Ri
                    TempLaunchAngle_inform = site4BM4 ;
                otherwise
                    disp('there are some problems!!(Launch Site)')        
            end
        otherwise
                disp('there are some problems!!(BM Type)')        
    end
        
    switch target(id).destinationid     
        case 1
            target(id).destination= struct('LatH','N','LatDeg',37.53,'LonH','E','LonDeg',126.96)    ; % seoul
            target(id).launchDirec= struct('AzDeg',TempLaunchAngle_inform(1,1)+0.3-0.6*rand(1,1),'ElDeg',TempLaunchAngle_inform(2,1)+0.004-0.008*rand(1,1))   ; 
        case 2 
            target(id).destination= struct('LatH','N','LatDeg',36.35,'LonH','E','LonDeg',127.38)    ; % daejeon
            target(id).launchDirec= struct('AzDeg',TempLaunchAngle_inform(1,2)+0.3-0.6*rand(1,1),'ElDeg',TempLaunchAngle_inform(2,2)+0.004-0.008*rand(1,1))   ; 
        case 3
            target(id).destination= struct('LatH','N','LatDeg',37.34,'LonH','E','LonDeg',127.92)    ; % wonju
            target(id).launchDirec= struct('AzDeg',TempLaunchAngle_inform(1,3)+0.3-0.6*rand(1,1),'ElDeg',TempLaunchAngle_inform(2,3)+0.004-0.008*rand(1,1))   ; 
        case 4
            target(id).destination= struct('LatH','N','LatDeg',35.87,'LonH','E','LonDeg',128.6)     ; % daegu
            target(id).launchDirec= struct('AzDeg',TempLaunchAngle_inform(1,4)+0.3-0.6*rand(1,1),'ElDeg',TempLaunchAngle_inform(2,4)+0.004-0.008*rand(1,1))   ; 
        case 5
            target(id).destination= struct('LatH','N','LatDeg',35.16,'LonH','E','LonDeg',126.9)     ; % gwangju
            target(id).launchDirec= struct('AzDeg',TempLaunchAngle_inform(1,5)+0.3-0.6*rand(1,1),'ElDeg',TempLaunchAngle_inform(2,5)+0.004-0.008*rand(1,1))   ; 
        case 6
            target(id).destination= struct('LatH','N','LatDeg',35.18,'LonH','E','LonDeg',129.08)    ; % busan
            target(id).launchDirec= struct('AzDeg',TempLaunchAngle_inform(1,6)+0.3-0.6*rand(1,1),'ElDeg',TempLaunchAngle_inform(2,6)+0.004-0.008*rand(1,1))   ; 
        otherwise
            disp('there are some problems!!(destination)')        
    end
     
end
