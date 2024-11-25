%%
%build theoretical mass list
charge=1.0072765*ESI_mode;
C13prob=0.011; %natural abundance of C13 
C12prob=1-0.011;
[e_array{1}, e_array{2}, e_array{3}, e_array{4}, e_array{5}, e_array{6}, e_array{7}, e_array{8}, e_array{9}, e_array{10}] = ndgrid(elements(1,1):elements(1,2),elements(2,1):elements(2,2),elements(3,1):elements(3,2),elements(4,1):elements(4,2),elements(5,1):elements(5,2),elements(6,1):elements(6,2),elements(7,1):elements(7,2),elements(8,1):elements(8,2),elements(9,1):elements(9,2),elements(10,1):elements(10,2));

el_matrix=[e_array{1}(:) e_array{2}(:) e_array{3}(:) e_array{4}(:) e_array{5}(:) e_array{6}(:) e_array{7}(:) e_array{8}(:) e_array{9}(:) e_array{10}(:)];
%% calculates m/z for theoretical assign table
formula_features=zeros(size(el_matrix,1),1);
valence_tot=zeros(size(el_matrix,1),1);
formula_features=formula_features+charge;
for els=1:no_els % all 10 elements
    formula_features(:,1)=formula_features(:,1)+el_matrix(:,els)*elements(els,3);%feature 1 is mass/charge and 3 is mass
    valence_tot=valence_tot+el_matrix(:,els)*elements(els,4);
end
valence_mod=mod(valence_tot,2); %valence number has to be even. zero if even 1 if odd. so only neutral 
%% calculates possibilities of elemental composition matrix
formula_features(:,2)=el_matrix(:,2)./(el_matrix(:,1)+el_matrix(:,10)); %feature 2 is H/C
formula_features(:,3)=el_matrix(:,3)./(el_matrix(:,1)+el_matrix(:,10)); %feature 3 is O/C
formula_features(:,4)=formula_features(:,1)*14/14.01565;   %Kendrick mass. reduction of mass
formula_features(:,5)=ceil(formula_features(:,4));  %Nom. Kendrick mass.round up value
formula_features(:,6)=formula_features(:,5)-formula_features(:,4);  %KMD
formula_features(:,7)=1+0.5*(2*el_matrix(:,1)+2*el_matrix(:,10)-el_matrix(:,2)-el_matrix(:,4)-el_matrix(:,6)); %DBE (1+ 0.5(2*C12+2*C13-H-N-P))
formula_features(:,8)=formula_features(:,7)-el_matrix(:,3);  %DBE-O
formula_features(:,9)=el_matrix(:,1)*C13prob/C12prob;  %predicted 13C pattern or 13C ratio
DBE_AIm=1+el_matrix(:,1)+el_matrix(:,10)-0.5*el_matrix(:,3)-el_matrix(:,5)-el_matrix(:,7)-0.5*(el_matrix(:,4)+el_matrix(:,2)+el_matrix(:,6));%DBE_AI 
C_AIm=el_matrix(:,1)+el_matrix(:,10)-0.5*el_matrix(:,3)-el_matrix(:,5)-el_matrix(:,4)-el_matrix(:,7)-el_matrix(:,6); %C_AI note el 7 is Se, why is this here??
formula_features(:,10)=DBE_AIm./C_AIm;  %Aromaticity index mod.
formula_features(formula_features(:,10)<0,10)=0;  %removes negative AImod
%formula_features(:,11)=formula_features(:,7)./(el_matrix(:,1)+el_matrix(:,10)); %feature 11 is DBE/C %created by MGD need to add a filter I think

formula=formula_features; %features 1-10 are characteristics of elemental compositions
elements_used=el_matrix;

% conditions: filters to remove bad formulas. %These definitions below define what to remove.
c1=formula(:,1)<min_mass | formula(:,1)>max_mass; %criteria 1 mass/charge filter. only accept mass/charge=150-800 
c2=formula(:,2)<0.3 | formula(:,2)>2.2; %remove if H/C<0.3 & H/C>2.2, accept only H/C=0.3-2.2
c3=formula(:,3)>1; %remove if O/C>1. accept only O/C= 0-1
c4=valence_mod==1; %remove if odd # of valence electrons (set equal to 1), i.e. only neutral electron valence accepted
c5=formula_features(:,8)>10|formula_features(:,8)<(-10); %DBE-O filter away if >10 or <-10 (Accept DBE=-10-10)
c6=el_matrix(:,4)>0&el_matrix(:,5)>0;% N and S cannot be together
c7=el_matrix(:,4)>0&el_matrix(:,10)>0; %N and C13 cannot be together
c8=el_matrix(:,5)>0&el_matrix(:,10)>0; %S and C13 cannot be together

conditions=c1+c2+c3+c4+c5+c6+c7+c8>0; %sums conditions and create big  filter

formula(conditions,:)=[]; %if condition above for c values is true (equal 1) then remove
elements_used(conditions,:)=[];

[formula,ind]=sortrows(formula,1);
elements_used=elements_used(ind,:);

refmass=formula(:,1);
formula_precut=formula;
elements_used_precut=elements_used;

clear C_AIm DBE_AIm valence_mod valence_tot
%% calibrate and assign
cd(myFolder) %locate folder with samples_all
sample_data=cell(nosamps,noreps);
for n=1:nosamps 
    n
    peakcandidates{n}=zeros(length(formula),noreps);
    assignment_errors{n}=zeros(length(formula),noreps);
    interference_errors{n}=zeros(length(formula),noreps);
    no_cand{n}=zeros(length(formula),noreps);
    for rep=1:noreps
        intensities_matrix{n,rep}=zeros(length(formula),1);
        clear mass_defect
        if data_type==1
            datsmpraw=table2array(readtable(names_samps{n,rep},'Delimiter',','));  %for csv files
        elseif data_type==2
            datx=mzxmlread(names_samps{n,rep});  %for mzXML files
            [peaks,time]=mzxml2peaks(datx,'Levels',1);   %for mzXML files
            datsmpraw=peaks{1};   %for mzXML files
        elseif data_type==3
            datsmpraw=xlsread(names_samps{n,rep});
        else
            cc=readtable(names_samps{n,rep});
            datsmpraw=[cc(:,2) cc(:,5) cc(:,3)];
            datsmpraw=table2array(datsmpraw);
        end

        raw_data{n,rep}=datsmpraw;
        full_length(n,rep)=length(datsmpraw);
        samplemax=max(datsmpraw(:,2)); %max intensity
        
        masses=datsmpraw(:,1);
        intensities=datsmpraw(:,2);
        mass_defects=datsmpraw(:,1)-floor(datsmpraw(:,1)); %value after decimal
        c13d2=(13.00335-12)/2;
        ppm_tol=2; %ppm tolerance
        flags=zeros(length(masses),1);
        %
        for m=1:length(masses)
            if mass_defects(m)>0.4&mass_defects(m)<0.8 % mass defect values 
                flags(m)=1; %masses within 0.4-0.8 get flag=1
                mono_mass=masses(m)-c13d2; %subtract C-13 mass difference/2
                m_diff=abs(masses-mono_mass); %subtract mass from new mass (w/o C-13)
                mono_pos=find(m_diff==min(m_diff)); %find  peak match if removed c-13
                if min(m_diff)/m*1e6<ppm_tol %max ppm
                    if intensities(mono_pos)/intensities(m)<10
                        flags(mono_pos)=2; %mass matched with mass defect peak get flag=2
                    end
                end
            end
        end
        data_filter=datsmpraw(flags==0,:);
        
        %noise removal
        noise{n,rep}=0;
        ken_mass=data_filter(:,1)*14/14.01565;
        ken_mass_defect=ceil(ken_mass)-ken_mass;
        %KMD slice method from MFAssign
        good_low=0.0011232*ken_mass+0.05;
        good_high=0.0011232*ken_mass+0.2;
        noise_mass=ken_mass_defect>good_low&ken_mass_defect<good_high;
        noise_peaks=data_filter(noise_mass,:);
        if length(noise_peaks)>100
            noise_quant=1*quantile(noise_peaks(:,2),0.99); %define noise
        else
            noise_quant=0;
        end
        noise{n,rep}=noise_quant;
 
        realpeaks=data_filter(:,2)>noise{n,rep};
        number_candidates(n,rep)=sum(realpeaks);
    
        datsmp=data_filter(realpeaks,:);
        %crude cal
            tmp_ppm=abs(data_filter(:,1)-cal_peak)./data_filter(:,1)*1e6;
            calposs=data_filter(tmp_ppm<linear_cal_ppm,:);
            calpk=calposs(calposs(:,2)==max(calposs(:,2)),1);
            ppmshift{n,rep}=(cal_peak-calpk)/cal_peak*1E6;
        
            if calchoice==1
                if abs(ppmshift{n,rep})>0
                datsmp(:,1)=datsmp(:,1)+ppmshift{n,rep}*datsmp(:,1)/1e6;
                else
                datsmp(:,1)=datsmp(:,1)+ppmshift{n-1,rep}*datsmp(:,1)/1e6;
                end
           else
                datsmp(:,1)=datsmp(:,1);
            end
        
        sample_data{n,rep}=datsmp;
    
        % fine calibration
        Kend=mod(formula(:,5),2)==1&elements_used(:,3)==formula(:,7)|elements_used(:,3)==formula(:,7)-1;
        noKend=Kend==0;
        neg_Kend=refmass;

        breakdown=round(length(datsmp)/100,0);
        d=zeros(length(datsmp),1);
        if breakdown==0
            breakdown=1;
        end
        for i=breakdown:length(datsmp)-breakdown
            if datsmp(i,2)>0.5*max(datsmp(i-(breakdown-1):i+(breakdown-1),2))%
                d(i)=1;
            end
        end
        
        datsmptest=datsmp(d==1,:);
        testint=zeros(length(formula),1);
        testerr=zeros(length(formula),1);
            for m=1:size(datsmptest,1)
                mass=datsmptest(m,1);
                intensity=datsmptest(m,2);
                tmp = abs(mass - neg_Kend);
                [~, idx] = min(tmp);
                ppmbrut=(neg_Kend(idx)-mass)/mass*1E6;
                ppmdiff=abs(ppmbrut);
                if ppmdiff<=fine_cal_ppm
                    testint(idx) = intensity;
                    testerr(idx) = ppmbrut;
                end
            end   
        testint(Kend==0)=0;
        testerr(Kend==0)=0;
    %
        calsmp{n,rep}(:,1)=formula(:,1);
        calsmp{n,rep}(:,2)=testint;
        calsmp{n,rep}(:,3)=testerr;
        if size(calsmp{n,rep},1)>0
            c=calsmp{n,rep}(:,2)>0.001*samplemax;
            cal_polynomial{n,rep}= polyfit(calsmp{n,rep}(c,1),calsmp{n,rep}(c,3),5);
        end
        %
        %assignment
        for m=1:size(sample_data{n,rep},1)
            mass1=sample_data{n,rep}(m,1);
            intensity=sample_data{n,rep}(m,2);
            
            % cal application
            if calchoice==1
                mass=mass1+polyval(cal_polynomial{n,rep},mass1)*mass1/1e6;
            else
                mass=mass1;
            end
            
            tmp = abs(mass - refmass);
                tmp_ppm=tmp/mass*1e6;
                
                mass_options=refmass(tmp_ppm<assign_ppm);
                option_diffs=abs(mass_options-mass);
                [~, idx] = min(tmp);
                
                ppmbrut=(refmass(idx)-mass)/mass*1E6;
                ppmdiff=abs(ppmbrut);
                if ppmdiff<=assign_ppm
                    intensities_matrix{n,rep}(idx) = intensity;
                    assignment_errors{n}(idx,rep) = ppmbrut;
                    no_cand{n}(idx,rep)=length(mass_options);
                end
         end
        
    end
end
%clear breakdown c calpk calposs charge d datsmp data_filter datsmptest
%clear full_length FWHM FWTM FWTMpeaks good_high good_low i idx idx_2nd
%clear intensity ken_mass ken_mass_defect Kend m mass mass1 n neg_Kend
%clear noise_mass noise_peaks noise_quant noKend ppmbrut ppmbrut_2nd ppmdiff
%clear realpeaks refmass rep samplemax testerr testint tmp tmp_ppm
%% Calibration graph
cd(HostFolder)
%shows calibration for first 3 samples, rep1
Figcal = figure('name','Calibration');
Figcal.Color=([1 1 1]);
Figcal.Units='centimeters';
Figcal.Position=([3 3 20 12]);
colz=2;
if nosamps>10
    plotno=10;
    rowz=5;
else
    plotno=nosamps;
    rowz=ceil(nosamps/2);
end

    for smp=1:plotno
        rep=1;
        subplot(rowz,colz,smp)
            c=calsmp{smp,rep}(:,2)>0;
            x1=linspace(150,800);
            y1=polyval(cal_polynomial{smp,rep},x1);
            hold on
            scatter(calsmp{smp,rep}(c,1),calsmp{smp,rep}(c,3),5,log10(calsmp{smp,rep}(c,2)),'filled')
            plot(x1,y1,'r','LineWidth',1)
        xlabel('m/z');ylabel('error (ppm)')
        legend('off')
        title(sampcodes{smp})
        axis([150 850 -2.5 2.5])
        calcb=colorbar();
    end
hold off
clear x1 y1 rep smp c
%% sample intensities

sample_intensities=cell(1,nosamps);
all_intensities=zeros(length(formula_precut),noreps*nosamps);

clear smp_common
sum_all=zeros(length(formula_precut),1);

for n=1:nosamps
    sample_intensities{1,n}=zeros(length(formula_precut),noreps);
    sum_samp{n}=zeros(length(formula_precut),1);
    for rep=1:noreps
        
        iter=(n-1)*noreps+rep;
        dat_out=intensities_matrix{n,rep};
        blnk=zeros(length(formula),1);
        %blnk=(intensities_matrix{69,1}+intensities_matrix{70,1})/2; %blank
        dat_out(dat_out<5*blnk)=0;
        sample_intensities{n}(:,rep)=dat_out;
        all_intensities(:,iter)=dat_out;
        peakfind=dat_out>0;
        sum_all=sum_all+peakfind;
        sum_samp{n}=sum_samp{n}+peakfind;
        
    end
end
sumx=sum(all_intensities(:,1:9)); %
%% Filters for peaks
refilter=sum_all>0 & elements_used(:,10)==0; %filters peaks to be in at least x samples (change from sum_all>0 to 1), no isotopologues
sum(refilter)
formula=formula_precut(refilter,:);
elements_used=elements_used_precut(refilter,:);
sum_all=sum_all(refilter,:);
all_intensities=all_intensities(refilter,:);
for n=1:nosamps
    sample_intensities{n}=sample_intensities{n}(refilter,:);
    peakcandidates{n}=peakcandidates{n}(refilter,:);
    assignment_errors{n}=assignment_errors{n}(refilter,:);
    interference_errors{n}=interference_errors{n}(refilter,:);
end
NOSC=4-(((4*(elements_used(:,1)+elements_used(:,10)))+(1*elements_used(:,2))-(3*elements_used(:,4))-(2*elements_used(:,3))-(2*elements_used(:,5))+(5*elements_used(:,6))+(2*elements_used(:,7))-1)./(elements_used(:,1)+elements_used(:,10)));
%% replicate, sum and filter
for smp=1:nosamps
    smp_detected{smp}=sample_intensities{smp}>0;
    smp_detected_no{smp}=sum(smp_detected{smp},2);
    smp_confident(:,smp)=smp_detected_no{smp}==noreps;
    for rep=1:noreps
        clean_rep=sample_intensities{smp}(:,rep).*smp_confident(:,smp);
        norm_rep(:,rep)=clean_rep/sum(clean_rep)*1e6;%normalize intensities by dividing each intensity in a sample by the sum of all non-normalized intensities for that sample them multiply to 1 million. this makes the sum of avg intensities equal to 1e6
    end
    avg_intensities(:,smp)=mean(norm_rep,2);
end
%
save('matlab') %create a file with all variables