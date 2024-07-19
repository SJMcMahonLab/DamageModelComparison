% Simple DNA model, reads in chemistry or physics-only model tuple data  
% and uses energy depositions in DNA strands, converting to SSBs by energy-dependent
% random sampling, clustering into DSBs, computing DSB complexity
clc;
close all;
seed = rng(1);

%Nucleus and DNA Parameters (dimensions in micrometers)
Rnuc = 4.65;
Rdna = 0.00115;
BPlength = 0.00034;
DNAbp = 6*10^9;

BPvolume = pi*(Rdna)^2*BPlength;
DNAvolume = DNAbp*BPvolume;
Nuclearvolume =4/3*pi*(Rnuc)^3;
DNAvolumefraction = DNAvolume/Nuclearvolume;
DNAarea = pi*(Rdna)^2;
Strandseparation = sqrt((DNAarea*4/3*pi*(Rnuc)^3)/DNAvolume);

FirstSeed = 1;
LastSeed = 50;

SSBs = 850;
DSBSeparation = 0.0034;
ComplexityLength = 0.0034;
Record_DamageSummary = [];

for Seed = FirstSeed:LastSeed
    disp(Seed)
    filename = sprintf('%s%d%s','NucleusTuple_HS_5MeV_Seed', Seed, '.phsp');
    pTuple = importdata(filename);
    
    if isempty(pTuple) 
        DamageSummary_Data = [NaN, NaN, NaN];
                
    elseif ~isempty(pTuple) 
    
        eventID = str2double(pTuple.textdata(:,5));
        xdata = str2double(pTuple.textdata(:, 2));
        ydata = str2double(pTuple.textdata(:, 3));
        zdata = str2double(pTuple.textdata(:, 4));
        energydepositions = pTuple.data(:,2); 
        dYList = zeros(size(energydepositions,1),1);
        t = [eventID, xdata, ydata, zdata, energydepositions,dYList];
        clearvars xdata ydata zdata energydeposition pTuple

        Record_DSB_Data = [];

        for interaction = 1:size(energydepositions,1)
            %Rotate grid by moving X,Y,Z 
            X = t(interaction,2);
            Y = t(interaction,3);
            Z = t(interaction,4);
            %distance to nearest strand
            dX = X - Strandseparation*round(X/Strandseparation);
            dY = Y - Strandseparation*round(Y/Strandseparation);
            dZ = Z - Strandseparation*round(Z/Strandseparation);
            theta = -(mod((floor(Z/BPlength)),10)*pi/5.0);
            dYPrime = dX*sin(theta) + dY*cos(theta);

            R = sqrt(dX^2+dY^2);

            %if outside DNA radius set energy to zero
                if R > Rdna
                    t(interaction,5) = 0;
                else
                    t(interaction,6) = dYPrime;
                end
        end
        %remove outside dna interactions 
        DNAinteractions = t(logical(t(:,5)),:);
        EnergyinDNA = sum(DNAinteractions(:,5));    %keV
        SSBEnergy = EnergyinDNA/SSBs;
        energyevents = numel(DNAinteractions(:,5));
        repeats = 1000; 
        RecordSSB = 0;
        RecordDSB = 0;
        RecordRepeat = [];
        OutputDamageData = [];

        for r = 1:repeats
            SSBCount = 0;
            DSBCount = 0;
            DSBPosition_List =[];
            HistoryID_List = [];
            RecordDamageData = [];
            SSBprob = zeros(energyevents,1);
            breakpositionx = zeros(energyevents,1);
            breakpositiony = zeros(energyevents,1);
            breakpositionz = zeros(energyevents,1);
            history = zeros(energyevents,1);
            Strand = zeros(energyevents,1);

            for event = 1:energyevents
                HistoryID = DNAinteractions(event,1);
                x = DNAinteractions(event,2);
                y = DNAinteractions(event,3);
                z = DNAinteractions(event,4);
                energy = DNAinteractions(event,5);
                SBprobability = energy/SSBEnergy;

                %convert into strand breaks
                if  SBprobability > rand 
                    SSBCount = SSBCount + 1;
                    breakpositionx(event,:) = x;
                    breakpositiony(event,:) = y;
                    breakpositionz(event,:) = z;
                    history(event,:) = HistoryID;
                    SSBprob(event,:) = SBprobability;   

                    %assign strand if dYPrime greater than zero             
                    if DNAinteractions(event,6) >= 0
                        Strand(event,:) = 1;
                    else
                        Strand(event,:) = 2;
                    end    
                end
            end

            SB_List = [breakpositionx, breakpositiony, breakpositionz, SSBprob, Strand, history];
            SB_List = SB_List(logical(SB_List(:,4)),:);

            TotalSB = size(SB_List,1);
            Strand1 = zeros(TotalSB,6);
            Strand2 = zeros(TotalSB,6);
            Break1 = [];
            Break2 = [];

            %separating SSB list into list for each strand
            for SSB = 1:TotalSB
                if  SB_List(SSB,5) == 1
                    Strand1(SSB,:) = SB_List(SSB,:);
                else
                    SB_List(SSB,5) = 2;
                    Strand2(SSB,:) = SB_List(SSB,:);
                end
            end

            Strand1 = Strand1(logical(Strand1(:,4)),:);
            Strand2 = Strand2(logical(Strand2(:,4)),:);
            Strand1_SBCount = size(Strand1,1);
            Strand2_SBCount = size(Strand2,1); 

            %clustering into DSBs
            for j = 1:Strand1_SBCount
                Break1 = Strand1(j,:);
                History1 = Break1(:,6);
                x1 = Break1(:,1);
                y1 = Break1(:,2);
                z1 = Break1(:,3);
                isDSB = false;

                for k = 1:Strand2_SBCount
                    Break2 = Strand2(k,:);
                    x2 = Break2(:,1);
                    y2 = Break2(:,2);
                    z2 = Break2(:,3);

                    dx = x1-x2;
                    dy = y1-y2;
                    dz = z1-z2;
                    Separation = sqrt(dx.*dx + dy.*dy + dz.*dz);

                    if Separation <= DSBSeparation
                        isDSB = true;
                        Strand1_Position = [x1, y1, z1];
                        Strand2_Position = [x2, y2, z2];
                        %will record position of last matching strand 2 break 
                        Strand2(k,1)= -999;
                        Strand1(j,1)= -999999;
                    end  
                end

                if isDSB == true
                    DSBCount = DSBCount+1;
                    DSBPosition = [(Strand1_Position + Strand2_Position)/2];
                    DSBPosition_List = [DSBPosition_List; DSBPosition];
                    HistoryID_List = [HistoryID_List; History1];
                    
                end 
            end

            DSBComplexityList = size(DSBPosition_List,1);
            DSBRepeatList = size(DSBPosition_List,1);
            
            %Loop through SB list for DSB complexity
            for DSB_Complex = 1:DSBCount
                DSB_Placement = DSBPosition_List(DSB_Complex,:);
                DSBComplexityCount = 0;
                
                for SB_Complex = 1:SSBCount
                    SB_Placement = SB_List(SB_Complex,1:3);
                    X_separation = DSB_Placement(1)-SB_Placement(1);
                    Y_separation = DSB_Placement(2)-SB_Placement(2);
                    Z_separation = DSB_Placement(3)-SB_Placement(3);
                    DamageSeparation = sqrt((X_separation)^2+(Y_separation)^2+(Z_separation)^2);
    
                    if DamageSeparation <= ComplexityLength
                        DSBComplexityCount = DSBComplexityCount+1;
                        SB_List(SB_Complex,1:3) = [999,999,DSB_Complex];
                                              
                    end
                end
                %Set 2 SBs minimum, some DSB SBs included in other DSB
                %complexity 
                if DSBComplexityCount >= 2
                    DSBComplexityList(DSB_Complex,:) = DSBComplexityCount;
                else
                    DSBComplexityList(DSB_Complex,:) = 2;
                end
                DSBRepeatList(DSB_Complex,:) = r;
            end
            
            %add strand breaks for each repeat
            RecordSSB = RecordSSB + SSBCount;
            RecordDSB = RecordDSB + DSBCount;  
            
            RecordDamageData = [DSBRepeatList,HistoryID_List,DSBPosition_List,DSBComplexityList];
            RecordDamageData = sortrows(RecordDamageData,2);
            OutputDamageData = [OutputDamageData;RecordDamageData];

        end     
        
        ExposureID = zeros(RecordDSB,1);
        ExposureID = [OutputDamageData(:,1:2),ExposureID];
        ExposureID(1,3) = 2;
        
        for dsb = 2:RecordDSB
            if ExposureID(dsb,2) ~= ExposureID(dsb-1,2)
                ExposureID(dsb,3) = 1;
            end 
            if ExposureID(dsb,1) ~= ExposureID(dsb-1,1)
                ExposureID(dsb,3) = 2;
            end

        end
                

        Damage_Output = table(ExposureID(:,3), OutputDamageData(:,3:6));
        filename_1 = sprintf('%s%d', 'SimpleDNA_DamageOutput_05MeV_Seed', Seed);
        writetable(Damage_Output, filename_1, 'WriteVariableNames',0,'Delimiter','\t');
        
        SSB_Yield = RecordSSB/repeats;
        DSB_Yield = RecordDSB/repeats;
        DSB_Complexity = sum(OutputDamageData(:,6))/repeats;

        DamageSummary_Data = [SSB_Yield, DSB_Yield, DSB_Complexity];
        
    end    
    
    Record_DamageSummary = [Record_DamageSummary; DamageSummary_Data];
end

Average_DamageYields = mean(Record_DamageSummary, 1, 'omitnan');

Total_Seeds = (length(Record_DamageSummary(~isnan( Record_DamageSummary))))/3;
SEM_DamageYields = std(Record_DamageSummary, 1, 1,'omitnan')/sqrt(Total_Seeds);

Output_DamageSummary = table(Average_DamageYields,SEM_DamageYields);
filename_2 = sprintf('%s','SimpleDNA_DamageSummary_05MeV');
writetable(Output_DamageSummary, filename_2, 'WriteVariableNames',0,'Delimiter','\t');

%exit;



