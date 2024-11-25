%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% csv assign version 7
%
% Originally written by Jeffrey Hawkes
%
% Feel free to distribute, but please acknowledge the programmers and curators!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AEN DATA 
%% Define & load files 
clear all
%reassign=1; % choose 0 to use old data, 1 to apply the same

myFolder='/Users/mariagdigernes/OneDrive - NTNU/MATLAB/csv_assign_maria-AeN_MANUSCRIPT/samples_all';
HostFolder='/Users/mariagdigernes/OneDrive - NTNU/MATLAB/csv_assign_maria-AeN_MANUSCRIPT';

cd(HostFolder) %command to go to folder with all files 
data_type=2;% 1=csv, 2=mzXML, 3=xlsx, 4=txt
ESI_mode=-1; %negative mode ?
[smpdat,names_samps] = xlsread('samplenamesAeNMS.xlsx');% creates two new variables smpdat (values only) and namessamps (names only).
nosamps=size(names_samps,1); %number of samples. calculates the number of samples on samplename file and stores on column 1
noreps=1;
sampnames=names_samps(:,1); % creates samp names as variable for long sample name (located on column 1 of samplenames file)
sampcodes=names_samps(:,2);  % creates samp codes as variable for short sample name (located on column 2 of samplenames file)

%sampcodes={'ESFA';'PLFA';'SRFA'}  ;
%
cal_peak=[];
cal_peak=369.11911; %remove if you do not need a linear internal cal
linear_cal_ppm=10;
fine_cal_ppm=3;
assign_ppm=0.7;
calchoice=0;

min_mass=150; %used this in assign.m limits for mass/charge=150-800
max_mass=800;

%elements: [min,max,exact mass,valence] this creates a table with each
%value specified below for each column
elements(1,:)=[4,50,12,4]; %C 4-50
elements(2,:)=[4,100,1.007825,1]; %H 4-100
elements(3,:)=[2,40,15.9949146,2]; %O 2-40
elements(4,:)=[0,2,14.003074,3]; %N 0-2
elements(5,:)=[0,1,31.972071,4]; %S 0-1
elements(6,:)=[0,0,30.973762,5]; %P
elements(7,:)=[0,0,79.916521,6]; %Se
elements(8,:)=[0,0,22.989770,1]; %Na
elements(9,:)=[0,0,34.96885269,7]; %Cl
elements(10,:)=[0,1,13.00335,4]; %C13

no_els=size(elements,1); %counts total number of elements and puts the number on column 1

assign

%if reassign==1
%assign %runs script called assign 
%else
%load('matlab') %runs all of the variables saved on matlab file (.mat)
%end

[smpdat,names_samps] = xlsread('samplenamesAeNMS.xlsx'); %

nosamps=size(names_samps,1);%number of samples. calculates the number of samples on samplename file and stores on column 1
noreps=1;
sampnames=names_samps(:,1);
sampcodes=names_samps(:,2);  
samplenamesfilter=9; %column to plot or current working column from samplenames file

%% PLFA
if false 
%PLFAmean and SD
PLFA_dist=triu(dist_matrix(116:121,116:121));
PLFA_mean_dist=mean(PLFA_dist(PLFA_dist>0))
PLFA_sd_dist=std(PLFA_dist(PLFA_dist>0))


PLFA_dist=triu(dist_matrix(116:121,116:121));
PLFA_mean_dist=mean(PLFA_dist(PLFA_dist>0))
PLFA_sd_dist=std(PLFA_dist(PLFA_dist>0))
%
PLFA=readtable("common_PLFA_neg.csv");
PLFA_mz=table2array(PLFA(:,7));
for x=1:length(formula)
   
PLFAmatch(x,1)=ismember(round(formula(x,1),4),round(PLFA_mz,4));
logi=PLFAmatch==1;   
end
PLcount=0;
for smpx=116:121
    PLcount=PLcount+1;
    PLFA_case=all_intensities(:,smpx);
    countPLFA=sum(PLFAmatch==1&PLFA_case>0);
    percPLFA=countPLFA/length(PLFA_mz)*100;
   
    PLFA_tab(PLcount,1)=countPLFA;
    PLFA_tab(PLcount,2)=percPLFA;
end
%
smp=0;
for smpx=116:121
    smp=smp+1;
    OC_PLFA(smp,1)=sum(all_intensities(logi,smp).*formula((logi),3))/sum(all_intensities(logi,smp));
    HC_PLFA(smp,1)=sum(all_intensities(logi,smp).*formula((logi),2))/sum(all_intensities(logi,smp));
    MW_PLFA(smp,1)=sum(all_intensities(logi,smp).*formula((logi),1))/sum(all_intensities(logi,smp));
end
 
table_PLFA=[mean(OC_PLFA) std(OC_PLFA); mean(HC_PLFA) std(HC_PLFA); mean(MW_PLFA) std(MW_PLFA)];
end 

%% Table 2 
if true
form_rnd=round(formula(:,8),4);


for smp=1:nosamps
    nopeaks(smp,1)=sum(avg_intensities(:,smp)>0); %col 1
    tot_int=sample_intensities{smp}(smp_confident(:,smp),:); %sample intensities
    TAC(smp,1)=sum(tot_int(:)); %col 2 TAC= sum of all_intensities per sample. originally divided by 2??
   
    OC_weighted_average=sum(avg_intensities(:,smp).*formula(:,3))/sum(avg_intensities(:,smp));
    OC_wa(smp,1)=OC_weighted_average;
    OC_weighted_variance=sum(avg_intensities(:,smp).*(formula(:,3)-OC_weighted_average).^2)/sum(avg_intensities(:,smp));
    OC_wvar(smp,1)=OC_weighted_variance;
    OC_weighted_std_deviation = sqrt(OC_weighted_variance);
    OC_wstdev(smp,1) = OC_weighted_std_deviation ;
    
    HC_weighted_average=sum(avg_intensities(:,smp).*formula(:,2))/sum(avg_intensities(:,smp));
    HC_wa(smp,1)=HC_weighted_average;
    HC_weighted_variance=sum(avg_intensities(:,smp).*(formula(:,2)-HC_weighted_average).^2)/sum(avg_intensities(:,smp));
    HC_wvar(smp,1)=HC_weighted_variance;
    HC_weighted_std_deviation = sqrt(HC_weighted_variance);
    HC_wstdev(smp,1) = HC_weighted_std_deviation ;

    MW_weighted_average=sum(avg_intensities(:,smp).*formula(:,1))/sum(avg_intensities(:,smp));%weighted average MW takes into account how high intensity a compound is. sum and divide by average. 
    MW_wa(smp,1)=MW_weighted_average;
    MW_weighted_variance=sum(avg_intensities(:,smp).*(formula(:,1)-MW_weighted_average).^2)/sum(avg_intensities(:,smp));
    MW_wvar(smp,1)=MW_weighted_variance;
    MW_weighted_std_deviation = sqrt(MW_weighted_variance);
    MW_wstdev(smp,1) = MW_weighted_std_deviation ;
    
    AImod_weighted_average=sum(avg_intensities(:,smp).*formula(:,10))/sum(avg_intensities(:,smp));
    AImod_wa(smp,1)=AImod_weighted_average;
    AImod_weighted_variance=sum(avg_intensities(:,smp).*(formula(:,10)-AImod_weighted_average).^2)/sum(avg_intensities(:,smp));
    AImod_wvar(smp,1)=AImod_weighted_variance;
    AImod_weighted_std_deviation = sqrt(AImod_weighted_variance);
    AImod_wstdev(smp,1) = AImod_weighted_std_deviation ;
    
    cut = 0; 
    datplot=avg_intensities(:,smp); 
    b=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==0;% N=0 S=0,  number of CHO peaks col 6
    CHO(smp,1)=sum(b);
    CHO_relintensity(smp,1) = sum(datplot(b))./sum(datplot); %gives the decimal for average relative abundance of CHO in each sample (i.e. take sum of CHO intensities and divide by the total abundance of CHO for each sample)
    
    c=datplot>cut&elements_used(:,4)>0&elements_used(:,5)==0; % requires N>0 S=0, number of CHON peaks col 7
    CHON(smp,1)=sum(c);
    CHON_relintensity(smp,1) = sum(datplot(c))./sum(datplot);
    
    d=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==1; %N=0 S=1, number of CHOS peaks 
    CHOS(smp,1)=sum(d);
    CHOS_relintensity(smp,1) = sum(datplot(d))./sum(datplot);
    
    DBE(smp,1)=sum(avg_intensities(:,smp).*formula(:,7))/sum(avg_intensities(:,smp));% DBE weighted average col 9 
    DBEminO(smp,1)= sum(avg_intensities(:,smp).*formula(:,8))/sum(avg_intensities(:,smp));% DBE-O weighted average col 10  

clear  b c d cut 
table2=[nopeaks TAC OC_wa OC_wstdev HC_wa HC_wstdev MW_wa AImod_wa AImod_wstdev CHO CHON CHOS DBE DBEminO CHO_relintensity CHON_relintensity CHOS_relintensity]; 
end
end


%% goodsmp badsmp
goodsmp=1:nosamps;
badsmp = smpdat(:,3)<10|smpdat(:,1)==0|table2(:,1)<800; % 
goodsmp(badsmp)=[]; %removes badsamples
sampcodes_cut=sampcodes(goodsmp);

%% Mass spectra
FigMS = figure('Name','Mass spectra'); %
FigMS.Color=([1 1 1]);
FigMS.Units='centimeters';
FigMS.Position=([10 10 30 30]);
%
sample_startval_ms = 1 ;%added start sample number 
for i=1:8  %position to graph
        rep=1;
        subplot(4,2,i) %(rows,columns,position) number of plots 
        smp= i + sample_startval_ms - 1 ;%changed added to connect sample start val to smp
        datplot=avg_intensities(:,smp);
        stem(formula(:,1),datplot,'Marker','none') %formula contains 
        axis([150 850 0 prctile(datplot(datplot>0),99.9)])
        ylabel('ion abundance'),xlabel('m/z')
        title(sampcodes{smp});      
end 

%% Fig.S3 overlay mass spectra
FigMSOverlay = figure('Name','Mass spectra overlay');%
FigMSOverlay.Color=([1 1 1]);
FigMSOverlay.Units='centimeters';
FigMSOverlay.Position=([10 5 30 10]);


smp_under1 = 59; %19=Q1P6822 HMW 59=Q4p1-150m
smp_over1 = 53; %66:Q4p410m. 45=Q2_P6_60m LMW
datplot_under1 = avg_intensities(:,smp_under1);
datplot_over1 = avg_intensities(:,smp_over1);

subplot(2,2,1); %number of plots 
stem(formula(:,1),datplot_under1,'Marker','none','Color',[0 0.4470 0.7410],'LineWidth',2);
hold on
stem(formula(:,1),datplot_over1,'Marker','none','Color',[0.9290 0.6940 0.1250],'LineWidth',2);
axis([150 850 0 prctile(datplot_under1(datplot_under1>0),99)]) %xmin xmax ymin ymax
ylabel('ion abundance'),xlabel('m/z')

legend([num2str(sampcodes{smp_under1})],[num2str(sampcodes{smp_over1})]);
hold on

subplot(2,2,2);
smp_under2 = 53; %LMW
smp_over2 = 59; %HMW  last one you plot comes on top
datplot_under2 = avg_intensities(:,smp_under2);
datplot_over2 = avg_intensities(:,smp_over2);
stem(formula(:,1),datplot_under2,'Marker','none','Color',[0.9290 0.6940 0.1250],'LineWidth',2);
hold on
stem(formula(:,1),datplot_over2,'Marker','none','Color',[0 0.4470 0.7410],'LineWidth',2);
axis([150 850 0 prctile(datplot_under2(datplot_under2>0),99)]) %xmin xmax ymin ymax
ylabel('ion abundance'),xlabel('m/z')

legend([num2str(sampcodes{smp_under2})],[num2str(sampcodes{smp_over2})]);
%next sample
smp_under3 = 19; % HMW 
smp_over3 = 35; %  LMW
datplot_under3 = avg_intensities(:,smp_under3);
datplot_over3 = avg_intensities(:,smp_over3);

subplot(2,2,3); %number of plots 
stem(formula(:,1),datplot_under3,'Marker','none','Color',[0 0.4470 0.7410],'LineWidth',2);
hold on
stem(formula(:,1),datplot_over3,'Marker','none','Color',[0.9290 0.6940 0.1250],'LineWidth',2);
axis([150 850 0 prctile(datplot_under3(datplot_under3>0),99)]) %xmin xmax ymin ymax
ylabel('ion abundance'),xlabel('m/z')

legend([num2str(sampcodes{smp_under3})],[num2str(sampcodes{smp_over3})]);
hold on
subplot(2,2,4);
smp_under4 = 35; %LMW
smp_over4 = 19; %HMW 19=Q1P6822 last one you plot comes on top
datplot_under4 = avg_intensities(:,smp_under4);
datplot_over4 = avg_intensities(:,smp_over4);
stem(formula(:,1),datplot_under4,'Marker','none','Color',[0.9290 0.6940 0.1250],'LineWidth',2);
hold on
stem(formula(:,1),datplot_over4,'Marker','none','Color',[0 0.4470 0.7410],'LineWidth',2);
axis([150 850 0 prctile(datplot_under4(datplot_under4>0),99)]) %xmin xmax ymin ymax
ylabel('ion abundance'),xlabel('m/z')

legend([num2str(sampcodes{smp_under4})],[num2str(sampcodes{smp_over4})]);
 



%% Fig.S5 VK CHON overlay combo
if true 
    FigVKCHON = figure('Name','VK:CHON in Q4 Q1 Q2');
    FigVKCHON.Color=([1 1 1]);
    FigVKCHON.Units='centimeters';
    FigVKCHON.Position=([5 3 13 10]);
    cut=0;
    sample_startval_vk = 1 ; 
    for i=1:89 %position on subplot. also total number of samples. must set 0 to badsmp samples as these are not removed her. no filter applied here.
        smp = i + sample_startval_vk -1; %added to connect smp to sample startval
        if (smpdat(i,1) ~= 0)% samples not equal to 0 on smpdat are removed.  col 7 are 3 cruises on 1 plot.
            %subplot(3,1,smpdat(smp,1)) %plot sample number and column on samplenames file
            datplot=avg_intensities(:,smp); %to obtain plots from avg_intensities (normalized)
            sizex=datplot/max(datplot)*1000; % 
            %
            a=datplot>cut; %all intensities from selected samples
            b=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==0;% N=0 S=0,  CHO
            c=datplot>cut&elements_used(:,4)>0&elements_used(:,5)==0; % N>0 S=0, CHON
            d=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==1; %N=0 S=1, CHOS

            hold on
            %scatter(formula(a,3),formula(a,2),12,'filled') %O/C,H/C,   number after Y represents size of each dot. each column represents a characteristic of elemental compositions (H/C, DBE,13C,KMD etc)
            %scatter(formula(b,3),formula(b,2),sizex(b),'k','filled') %CHO 'r'red 'k' black etc. 
            %scatter(formula(c,3),formula(c,2),sizex(c),'r','filled') %CHON 
            %scatter(formula(d,3),formula(d,2),sizex(d),'g','filled') %CHOS
           
            if smpdat(smp,1) == 1 %Q1 
                %subplot(1,1,1)
                %scatter(formula(a,3),formula(a,2),sizex(a),"MarkerEdgeColor","#556B2F",,'LineWidth',3) %O/C,H/C,green color
                %scatter(formula(b,3),formula(b,2),sizex(b),"MarkerEdgeColor","#556B2F",'LineWidth',2) %CHO
                scatter(formula(c,3),formula(c,2),sizex(c),"MarkerFaceColor","#FF6666","MarkerEdgeColor","#FF6666") %CHON red
                %scatter(formula(d,3),formula(d,2),sizex(d),"MarkerEdgeColor","#556B2F",'LineWidth',3) %CHOS
            elseif smpdat(smp,1) == 2 %Q2
                %subplot(1,1,1)
                %scatter(formula(a,3),formula(a,2),sizex(a),"MarkerEdgeColor","#FFC20A",'LineWidth',2) %O/C,H/C,yellow gold color
                %scatter(formula(b,3),formula(b,2),sizex(b),"MarkerEdgeColor","#FFC20A",'LineWidth',2) %CHO
                scatter(formula(c,3),formula(c,2),sizex(c),"MarkerFaceColor","#8B4513", "MarkerEdgeColor","#8B4513") %CHON 
                %scatter(formula(d,3),formula(d,2),sizex(d),"MarkerEdgeColor","#FFC20A",'LineWidth',2) %CHOS
            elseif smpdat(smp,1)== 3 %Q4
                %subplot(1,1,1) 
                %scatter(formula(a,3),formula(a,2),sizex(a),"MarkerEdgeColor","#001f3F",'LineWidth',3)% dark blue color
                %scatter(formula(b,3),formula(b,2),sizex(b),"MarkerEdgeColor","#001f3F",'LineWidth',2) %CHO
                scatter(formula(c,3),formula(c,2),sizex(c),"MarkerFaceColor","#000000",'MarkerEdgeColor',"#000000") %CHON
                %scatter(formula(d,3),formula(d,2),sizex(d),"MarkerEdgeColor","#001f3F",'LineWidth',3) %CHOS
            end
            axis([0 1 0 2.4])
            ylabel('H/C'),xlabel('O/C')
            ax = gca;               % Get the current axis
            ax.FontSize = 14;       % Set font size for tick labels on both x-axis and y-axis
            %title(sampcodes{'t0', 't1'});
            %title(sampcodes{smp}); 
        end
    end
clear b c d cut 

colorbar('off');
colormap('default');
colorbars = findobj(FigVKCHON, 'Type', 'Colorbar');
if ~isempty(colorbars)
    delete(colorbars);  % Remove any existing colorbar
end

end

%% TAC TIC signals
if false
for i=1:nosamps
    TICx(i,1)=sum(sample_data{i,1}(:,2));
    TICn(i,1)=sum(raw_data{i,1}(:,2));
end
figure('Name','TAC TIC signals');
subplot(2,1,1)
hold on
plot(1:nosamps,TAC)
plot(1:nosamps,TICx)
plot(1:nosamps,TICn)
legend('TAC','signal','TIC')
subplot(2,1,2)
hold on
plot(1:nosamps,TAC./TICx)
legend('TAC/signal')
end
%% Bray curtis dissimilarity 
dist_matrix=zeros(nosamps,nosamps); %creates dissimilarity matrix template based on number of samples 
for dist_x=1:nosamps
    for dist_y=1:nosamps
        set1=avg_intensities(:,dist_x);
        set2=avg_intensities(:,dist_y);
        T1=set1/sum(set1); % T1 T2 are normalized here
        T2=set2/sum(set2);
        %T1=set1; %OBS! horseshoe effect occurs when not normalized T1 T2
        %T2=set2;
        T3=abs(T1-T2); %
        T4=T1+T2;
        diss=sum(T3(:))/sum(T4(:))*100; %B-C Diss formula used
        dist_matrix(dist_x,dist_y)=diss;
    end
end
clear T1 T2 T3 T4 diss


dist_avg=mean(dist_matrix,1);
mean_samp_no=find(dist_avg==min(dist_avg));
dist_matrix_cut=dist_matrix(goodsmp,goodsmp);
BCD_U=triu(dist_matrix_cut);
b=find(BCD_U~=0);
median(BCD_U(b))
quantile(BCD_U(b),0.75)
quantile(BCD_U(b),0.25)
max(BCD_U(b))

[Y,eigvals] = cmdscale(dist_matrix_cut);
hold off 

%% Heatmap
%made from  Bray curtis dissimilarity (BCD) eig values

hcb2=colorbar();
Figheatmap= figure('Name','Heatmap BCD');
xvalues = string(sampcodes(goodsmp));
h=heatmap(xvalues, xvalues, dist_matrix_cut);

%
%FigheatmapWM = figure()
%h= heatmap(smpdat(goodsmp,6), smpdat(goodsmp,19), 'Water mass','H/C','ColorVariable', 'TemperatureF');


%% Dendrogram
 
FigDend = figure('Name','Dendrogram'); 
FigDend.Color=([1 1 1]);
FigDend.Units='centimeters';
FigDend.Position=([5 5 12 12]);
colormap winter
[rows,columns] = size(dist_matrix_cut);
v_mat = [];
    for i = 1:rows-1
         v_mat = [v_mat dist_matrix_cut(i,i+1:columns)];
    end
    Z=linkage(v_mat);
    Clu = cluster(Z,'Cutoff',15,'Criterion','distance');
[H,T,outperm]=dendrogram(Z,0,'ColorThreshold',15);
    axis([0.5 length(goodsmp)+0.5 0 50])
    xticklabels(sampcodes_cut(outperm))
    xtickangle(90)
    ticklabels = get(gca,'XTickLabel');
 lineColours = cell2mat(get(H,'Color'));
 colourList = unique(lineColours, 'rows');

myColours = [230,186,0;
             127,127,127;
             10,35,140;
             176,0,0;
             158,182,72;
             79,129,189;
             209,99,9;
             4,122,146]/255;

%// Replace each colour (colour by colour). Start from 2 because the first colour are the "unclustered" black lines             
for colour = 2:size(colourList,1)
    %// Find which lines match this colour
    idx = ismember(lineColours, colourList(colour,:), 'rows');
    %// Replace the colour for those lines
    lineColours(idx, :) = repmat(myColours(colour-1,:),sum(idx),1);
end
%// Apply the new colours to the chart's line objects (line by line)
for linex = 1:size(H,1)
    set(H(linex), 'Color', lineColours(linex,:));
end  

pntCols=myColours;
  pntCols(1,:)=[0 0 0];
  pntCols(6,:)=myColours(4,:);
  pntCols(4,:)=myColours(2,:);
  pntCols(5,:)=myColours(1,:);
  pntCols(2,:)=myColours(5,:);
%% PCoA
if true 
FigPCOA = figure('Name','PCoA'); 
FigPCOA.Color = ([1 1 1]);
FigPCOA.Units = 'centimeters';
FigPCOA.Position = ([5 5 12 12]);

%colormap jet
%colormap(cmocean('thermal'));%cmocean is an oceonography color matrix. made colormap into a function by putting paranthesis after it
%colormap(cmocean('algae'));
%colormap(cmocean('dense'));
%colormap(cmocean('haline'));
%colormap(crameri('-bamako'))
colormap(crameri('-batlow'))

xcoordinatevalues= Y(:,1);
ycoordinatevalues= Y(:,2);

scatter(xcoordinatevalues,ycoordinatevalues,100,smpdat(goodsmp,19), 'filled');%changed to smpdat from number of peaks to include colors on number
%scatter(Y(:,1),Y(:,2),150, nopeaks(goodsmp),'filled'); % Y(:,2) represents 2nd coordinate on PCOA. 150 number is size of dot
%scatter(Y(:,1),Y(:,2),20,pntCols(Clu,:),'filled');
 
text((xcoordinatevalues)+0.2,(ycoordinatevalues)+0.01,sampcodes(goodsmp)); %text labels

xlabel(['PCoA 1: ' num2str(100*eigvals(1)/sum(eigvals),'%-10.0f') '%']);
ylabel(['PCoA 2: ' num2str(100*eigvals(2)/sum(eigvals),'%-10.0f') '%']); 



 cb=colorbar('northoutside'); %position of axis
 %c.Label.String = 'Bacterial Production ug C /Ld'; 
 c.Label.String = 'H/Cwa'; 
 %c.Label.String = 'Temperature (°C)'; 
 %c.Label.String = 'DOC (μM)'; 

 %caxis([1e6 5e6])

axis([min(xcoordinatevalues)*1.1 max(xcoordinatevalues)*1.1 min(ycoordinatevalues)*1.1 max(ycoordinatevalues)*1.1])


set(gca, 'FontSize', 12);
end

%% Pearson correlation 
clear cors
cors(:,1)=Y(:,1);
cors(:,2)=Y(:,2);


for corchoice=1:2 %change to 3
        for i=1:length(formula)
            [r,p]=corr(avg_intensities(i,goodsmp)',cors(:,corchoice));
            cor_peak(i,corchoice)=r;
            plistpeak(i,corchoice)=p;
            if p<0.001 %originally p<0.015
                cor_peak_sig(i,corchoice)=r;
            else
                cor_peak_sig(i,corchoice)=0;
            end
        end
    cor_peak(isnan(cor_peak))=0;
end

FigPearsonPCOA= figure('Name','Pearson correlation'); 
FigPearsonPCOA.Color=([1 1 1]);
FigPearsonPCOA.Units='centimeters';
FigPearsonPCOA.Position=([5 1 20 15]);

sz=mean(all_intensities,2); % 
b=sz>0;
c1=sz>0&abs(cor_peak_sig(:,1))>0; 
c2=sz>0&abs(cor_peak_sig(:,2))>0; %change cor peak to sig 2 or 3 
%
subplot(2,2,1)
scatter(formula(b,3),formula(b,2),sz(b)/max(sz)*50,cor_peak(b,1),'filled')
hcb=colorbar();
hcb.Title.String = "Pearson's PCoA1";
    axis([0 1 0 2])
    ylabel('H/C'); xlabel('O/C')
%
subplot(2,2,2)
scatter(formula(b,3),formula(b,2),sz(b)/max(sz)*50,cor_peak(b,2),'filled') 
hcb=colorbar();
%hcb.Title.String = "DOC (μM)";
%hcb.Title.String = "Temperature";
%hcb.Title.String = "TAC";
%hcb.Title.String = "Bacterial production";
hcb.Title.String = "Pearson's PCoA2";
    axis([0 1 0 2])
    ylabel('H/C'); xlabel('O/C')    
%
subplot(2,2,3)
scatter(formula(c1,3),formula(c1,2),sz(c1)/max(sz)*50,cor_peak_sig(c1,1),'filled')
hcb=colorbar();
hcb.Title.String = "Sig. Pearson's PCoA1";
    axis([0 1 0 2])
    ylabel('H/C'); xlabel('O/C')
%
subplot(2,2,4)
scatter(formula(c2,3),formula(c2,2),sz(c2)/max(sz)*50,cor_peak_sig(c2,2),'filled') % change cor_peak_sig second position.
hcb=colorbar();
%hcb.Title.String = "Sig. Pearson's DOC (μM)";
%hcb.Title.String = "Sig. Pearson's Temperature";
%hcb.Title.String = "Sig. number of peaks";
%hcb.Title.String = "Sig. Pearson's Chl a";
hcb.Title.String = "Sig. Pearson's PCoA2";

axis([0 1 0 2])
ylabel('H/C'); xlabel('O/C')       


%colormap(cmocean('matter'));
%colormap(crameri('-imola'))
colormap(crameri('bamako'))



%% Fig.5 Boxplot OC

groupingvariable1 = smpdat(goodsmp,7); %7 is water masses
font_size_labels=28;
font_size_axis=36;


if true
    %Boxplot for OC and water masses all 3 cruises combined
    FigBoxplotPeak_OC = figure('Name','O/C and water masses', 'Position', [100, 100, 1500, 600]);% [left, bottom, width, height]
     % Set default font properties for the figure
    set(FigBoxplotPeak_OC, 'DefaultAxesFontName', 'Times New Roman');
    boxplot(smpdat(goodsmp,18),groupingvariable1,'Labels',{'AW/mAW','wPW','PW','AW/mAW','wPW','PW','EBDW','AW/mAW','wPW','PW','EBDW'});
    ax = gca;  % Get the current axis
    ax.YAxis.FontSize = font_size_labels;  % Increase font size for y-axis numbers (adjust as needed)
    % Set font size for x-axis labels (tick labels)
    ax.XAxis.FontSize = font_size_labels;  % Adjust the x-axis tick label font size
     % Increase font size for x-axis category labels by finding all text objects on the x-axis
    h = findobj(gca, 'Type', 'text');
    set(h, 'FontSize', font_size_axis);  % Set font size for each text object on the x-axis
    ylabel('w.a. O/C ratio', 'FontSize', font_size_axis)
    xlabel('early winter                         late winter                          spring  ', 'FontSize', font_size_axis)
     % Increase line thickness for all boxplot lines
    set(findobj(gca, 'Type', 'line'), 'LineWidth', 1.8);  % Adjust the line thickness (e.g., 1.5)
    ax.LineWidth = 1.8;  % Set the thickness of the plot's surrounding box
    ylim([0.440 0.489])
   
end 

if false
    %Boxplot for MW and water masses all 3 cruises combined
    FigBoxplotPeak_MW = figure('Name','MW and water masses', 'Position', [100, 100, 1500, 600]);
    set(FigBoxplotPeak_MW, 'DefaultAxesFontName', 'Times New Roman');
    boxplot(smpdat(goodsmp,21),groupingvariable1,'Labels',{'AW/mAW','wPW','PW','AW/mAW','wPW','PW','EBDW','AW/mAW','wPW','PW','EBDW'});
    ax = gca;  % Get the current axis
    ax.YAxis.FontSize = font_size_labels;  % Increase font size for y-axis numbers (adjust as needed)
    % Set font size for x-axis labels (tick labels)
    ax.XAxis.FontSize = font_size_labels;  % Adjust the x-axis tick label font size
     % Increase font size for x-axis category labels by finding all text objects on the x-axis
    h = findobj(gca, 'Type', 'text');
    set(h, 'FontSize', font_size_axis);  % Set font size for each text object on the x-axis
    ylabel('w.a. Molecular weight', 'FontSize', font_size_axis)
    xlabel('early winter                   late winter                        spring  ', 'FontSize', font_size_axis)
    % Increase line thickness for all boxplot lines
    set(findobj(gca, 'Type', 'line'), 'LineWidth', 1.8);  % Adjust the line thickness (e.g., 1.5)
    ax.LineWidth = 1.8;  % Set the thickness of the plot's surrounding box
    ylim([350 398])
   
end 


%% T-test
data = nopeaks(goodsmp);
groupA = smpdat(goodsmp,1)==3; %Q1 & Q2
groupB = smpdat(goodsmp,1)==1;%Q4
dataAHC = data(groupA);
dataBHC = data(groupB);
[h,p,ci,stats] = ttest2(dataAHC, dataBHC)

clear data groupA groupB dataA dataB
%ttest on intensities
if false
    data = Y(:,goodsmp);
    groupA = smpdat(goodsmp,1)==1;
    groupB = smpdat(goodsmp,1)==2;
    dataA = data(groupA);
    dataB = data(groupB);
    [h2,p2] = ttest2(dataA, dataB)
end 
%% Mann-Whitney-Wilcoxon Test H/C
total_data = smpdat(goodsmp,19); %H/C
sample_set1 = smpdat(goodsmp,1)==1;
sample_set2 = smpdat(goodsmp,1)==3;
grp1 = total_data(sample_set1);
grp2 = total_data(sample_set2);
X1=transpose(grp1);
X2=transpose(grp2);
mwwtest(X1,X2)


%% formula presence
% presence to absence of formulas (total 4652 formulas)

mzthreshold = 0;
Ngoodsamples = length(goodsmp);
avg_intensities_logical = avg_intensities(:,goodsmp) > mzthreshold;
sumdata = sum(avg_intensities_logical,2); % 0 if present in no samples, Ngoodsamples if present in all samples, etc..
presentInAll = (sumdata == Ngoodsamples);
NpresentInAll = sum(presentInAll); %number of formulas present in all samples
index = transpose(1:size(presentInAll,1));
presentInAllIndex = nonzeros(presentInAll.*index);


figure ('name','formulas:presence')
plot(sumdata);
figure ()
plot(sort(sumdata));

%% correlation graphs
if false

FigQ1MW = figure('Name','Q1');
xq1axis= smpdat(goodsmp(1:26),21); %MW
yq1axis= smpdat(goodsmp(1:26),19); %H/C
scatter(xq1axis, yq1axis, 'MarkerFaceColor', 'b');
ylabel('H/C wa');
xlabel('Molecular weight');
title('Late winter');
% Fit a line using linear regression
p = polyfit(xq1axis, yq1axis, 1);  % p(1) is the slope, p(2) is the intercept
yfit = polyval(p, xq1axis);  % Compute the fitted y-values
% Plot the correlation line
hold on;
plot(xq1axis, yfit, '-r', 'LineWidth', 1.5);  % Add the line of best fit in red
hold off;
% Calculate R-squared
yresid = yq1axis - yfit;  % Residuals
SSresid = sum(yresid.^2);  % Sum of squared residuals
SStotal = (length(yq1axis) - 1) * var(yq1axis);  % Total sum of squares
Rsq = 1 - SSresid / SStotal;  % R-squared value
% Display R-squared value on the plot
text(min(xq1axis), max(yq1axis), ['R^2 = ' num2str(Rsq)], 'FontSize', 18, 'Color', 'r');
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Axis font size and tick mark thickness




FigQ2MW = figure('Name','Q2');
xaxis= smpdat(goodsmp(27:54),21); %MW
y1axis= smpdat(goodsmp(27:54),38); %H/C
scatter(xaxis, y1axis, 'MarkerFaceColor', 'b');
ylabel('H/C wa');
xlabel('Molecular weight');
title('Spring');
% Fit a line using linear regression
p = polyfit(xaxis, y1axis, 1);  % p(1) is the slope, p(2) is the intercept
yfit = polyval(p, xaxis);  % Compute the fitted y-values
% Plot the correlation line
hold on;
plot(xaxis, yfit, '-r', 'LineWidth', 1.5);  % Add the line of best fit in red
hold off;
% Calculate R-squared
yresid = y1axis - yfit;  % Residuals
SSresid = sum(yresid.^2);  % Sum of squared residuals
SStotal = (length(y1axis) - 1) * var(y1axis);  % Total sum of squares
Rsq = 1 - SSresid / SStotal;  % R-squared value
% Display R-squared value on the plot
text(min(xaxis), max(y1axis), ['R^2 = ' num2str(Rsq)], 'FontSize', 18, 'Color', 'r');
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Axis font size and tick mark thickness



FigQ4MW = figure('Name','Q4');
xq4axis= smpdat([55:62,64,66:76,78:81,84,86:88],21); %MW
yq4axis= smpdat([55:62,64,66:76,78:81,84,86:88],19); %H/C
%xq4axis= smpdat(goodsmp(55:89),21); %MW
%yq4axis= smpdat(goodsmp(55:89),19); %H/C
scatter(xq4axis, yq4axis, 'MarkerFaceColor', 'b');
ylabel('H/C wa');
xlabel('Molecular weight');
title('Early winter');
% Fit a line using linear regression
p = polyfit(xq4axis, yq4axis, 1);  % p(1) is the slope, p(2) is the intercept
yfit = polyval(p, xq4axis);  % Compute the fitted y-values
% Plot the correlation line
hold on;
plot(xq4axis, yfit, '-r', 'LineWidth', 1.5);  % Add the line of best fit in red
hold off;
% Calculate R-squared
yresid = yq4axis - yfit;  % Residuals
SSresid = sum(yresid.^2);  % Sum of squared residuals
SStotal = (length(yq4axis) - 1) * var(yq4axis);  % Total sum of squares
Rsq = 1 - SSresid / SStotal;  % R-squared value
% Display R-squared value on the plot
text(min(xq4axis), max(yq4axis), ['R^2 = ' num2str(Rsq)], 'FontSize', 18, 'Color', 'r');
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Axis font size and tick mark thickness
end

%%
if false
Figpeaknumber = figure('Name', 'Peak Number');
xaxis = 1:89; % Original x-axis data
yaxispeak = smpdat([1:89], 20); % Generate peak numbers (assumed 89 values)

% Sort the y-axis values (peak numbers) and get indices
[yaxispeak_sorted, sortIdx] = sort(yaxispeak);

% Apply the sorting to the x-axis
xaxis_sorted = sortIdx; % Use sorting indices directly as x-axis

% Plot the sorted data
scatter(1:length(xaxis),yaxispeak_sorted,100, 'filled');
%scatter(xaxis_sorted, yaxispeak_sorted, 'MarkerFaceColor', 'b');
%xticks(1:length(xaxis));
ylabel('Peak Number');
xlabel('Sample (Sorted by Peak Number)');
title('Peak Number per Sample (Sorted)');
end

%% Assigned data file
form_head={'C','H','O','N','S','P','Se','Na','Cl','13C','m/z','H/C','O/C','KM','NKM','KMD','DBE','DBE-O','13Cratio','AImod'};
for smp=1:nosamps
    headx{1,smp}=strcat(sampcodes{smp});
end
    heads=[form_head headx];
    assign_dat=[elements_used formula avg_intensities];
    csvwrite_with_headers('assigned_data.csv',assign_dat,heads)

