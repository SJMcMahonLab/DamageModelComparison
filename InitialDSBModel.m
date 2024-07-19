% Initial DSB damage model, reads in chemistry or physics-only model tuple data  
% and uses energy depositions in DNA strands, converting to SSBs by energy-dependent
% random sampling, clustering into DSBs, computing DSB complexity
clc;
close all;
seed = rng(1);

FirstSeed = 1;
LastSeed = 50;

DSBperGy = 111;
TargetRadius = 4.65; % in um
DSBEnergy = (4/3*pi()*TargetRadius^3*6.242)/DSBperGy; %in keV, where 1 Gy = 6.242 keV/um^3
disp(DSBEnergy) 
DSBComplexity = 0.43; %fraction of complex to simple breaks
Record_DamageSummary = [];

for Seed = FirstSeed:LastSeed
    disp(Seed)
    filename = sprintf('%s%d%s','NucleusTuple_HS_05MeV_Seed', Seed, '.phsp');
    pTuple = importdata(filename);
    
    if isempty(pTuple) 
        DamageSummary_Data = [NaN, NaN];
                
    elseif ~isempty(pTuple) 
    
        eventID = str2double(pTuple.textdata(:,5));
        xdata = str2double(pTuple.textdata(:, 2));
        ydata = str2double(pTuple.textdata(:, 3));
        zdata = str2double(pTuple.textdata(:, 4));
        energydepositions = pTuple.data(:,2); 
        DNAinteractions = [eventID, xdata, ydata, zdata, energydepositions];
        clearvars xdata ydata zdata energydeposition pTuple

        energyevents = numel(DNAinteractions(:,5));
        repeats = 1000; 
            
        RecordDSB = 0;
        OutputDamageData = [];

        for r = 1:repeats
            DSBCount = 0;
            DSB_List = [];
            DSBComplexityList = [];
            DSBRepeatList = [];
            RecordDamageData = [];
            
            DSBprob = zeros(energyevents,1);
            breakpositionx = zeros(energyevents,1);
            breakpositiony = zeros(energyevents,1);
            breakpositionz = zeros(energyevents,1);
            history = zeros(energyevents,1);

            for event = 1:energyevents
                HistoryID = DNAinteractions(event,1);
                x = DNAinteractions(event,2);
                y = DNAinteractions(event,3);
                z = DNAinteractions(event,4);
                energy = DNAinteractions(event,5);
                DSBprobability = energy/DSBEnergy;

                %convert energy into DSBs
                if  DSBprobability > rand 
                    DSBCount = DSBCount + 1;
                    breakpositionx(event,:) = x;
                    breakpositiony(event,:) = y;
                    breakpositionz(event,:) = z;
                    history(event,:) = HistoryID;
                    DSBprob(event,:) = DSBprobability;    
                end
            end
            
            DSB_List = [breakpositionx, breakpositiony, breakpositionz, DSBprob, history];
            DSB_List = DSB_List(logical(DSB_List(:,4)),:);

            TotalDSBs = size(DSB_List,1);
            
            DSBComplexityList = zeros(TotalDSBs,1);
            DSBRepeatList = zeros(TotalDSBs,1);
           
            %Add DSB complexity 
            for DSBComplex = 1:TotalDSBs
                if rand > DSBComplexity
                    DSBComplexityList(DSBComplex,1) = 2;
                    DSBRepeatList(DSBComplex,1) = r;
                else 
                    DSBComplexityList(DSBComplex,1) = 3;
                    DSBRepeatList(DSBComplex,1) = r;
                end
            end
             
            %Add DSB information for each repeat
            DSB_List = [DSBRepeatList, DSB_List,DSBComplexityList];            
            RecordDSB = RecordDSB + DSBCount;  
            
            RecordDamageData = sortrows(DSB_List,6);
            OutputDamageData = [OutputDamageData;RecordDamageData];

        end     
        
        ExposureID = zeros(RecordDSB,1);
        ExposureID = [OutputDamageData(:,1),OutputDamageData(:,6),ExposureID];
        ExposureID(1,3) = 2;
        
        for dsb = 2:RecordDSB
            if ExposureID(dsb,2) ~= ExposureID(dsb-1,2)
                ExposureID(dsb,3) = 1;
            end 
            if ExposureID(dsb,1) ~= ExposureID(dsb-1,1)
                ExposureID(dsb,3) = 2;
            end
        end
                
        Damage_Output = table(ExposureID(:,3), OutputDamageData(:,2:4),OutputDamageData(:,7));
        filename_1 = sprintf('%s%d', 'InitialDSBDamage_DamageOutput_05MeV_Seed', Seed);
        writetable(Damage_Output, filename_1, 'WriteVariableNames',0,'Delimiter','\t');
        
        DSB_Yield = RecordDSB/repeats;
        DSB_Complexity = sum(OutputDamageData(:,7))/repeats;
        DamageSummary_Data = [DSB_Yield, DSB_Complexity];
        
    end    
    
    Record_DamageSummary = [Record_DamageSummary; DamageSummary_Data];
end

Average_DamageYields = mean(Record_DamageSummary, 1, 'omitnan');

Total_Seeds = (length(Record_DamageSummary(~isnan( Record_DamageSummary))))/3;
SEM_DamageYields = std(Record_DamageSummary, 1, 1,'omitnan')/sqrt(Total_Seeds);

Output_DamageSummary = table(Average_DamageYields,SEM_DamageYields);
filename_2 = sprintf('%s','InitialDSBDamage_DamageSummary_05MeV');
writetable(Output_DamageSummary, filename_2, 'WriteVariableNames',0,'Delimiter','\t');

%exit;



