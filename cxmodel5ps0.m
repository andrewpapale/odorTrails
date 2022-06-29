%
% Version of cxmodel with disordered parameters and activity transmission
% delays.
%
% nReqAvs is the number of requested avalanches per node
% Expect at most nReqAvs*N, depending on MaxDuration
%
% Optimizes data precision
%
% Rashid V. Williams-Garcia 12/2/15

function [cShapes,cWebs,rngState] = cxmodel5ps0(weightmat,Ncw,Taurs,sig_del,delSDs,MaxSteps)
%     rng('shuffle','twister');
    temp = rng;
    rngState = temp.Seed;
%     fname = strcat('cxm',num2str(rngState));
    
    N = length(weightmat);

    connmat = cell(N,1);
    for i=1:N
        connmat{i} = find(weightmat(i,:));
    end
    connects = cellfun(@numel,connmat);

    cShapes = cell(Ncw,1);
    cWebs = cell(Ncw,1);
    
    %% Begin Simulation
    if Ncw>1
        parpool
        parfor n=1:Ncw
            prevZ = zeros(N,1);%zeros(N,1,pTaur);
            t = 1;
            dt = 0; %#ok<NASGU>
            
            DrivenSites = [randi(N),t];   %[driven site,activation time]
            actives = DrivenSites(DrivenSites(:,2)==t);
            prevZ(actives) = 1;
            activations = DrivenSites;
            cSteps = 1;
            cWeb = zeros(0,2);
            pPairs =  zeros(0,3);%potential causal pairs
            %1st col gives presynaptic activation
            %2nd col gives postsynaptic site
            %3rd col gives postsynaptic activation time
            
            while ~isempty(DrivenSites) && cSteps<=MaxSteps
                %acts = [91 1;92 2;93 3; 94 3; 95 3];
                
                %determine activation indeces:
                [~,preActs] = ismember([actives t*ones(size(actives))],activations,'rows');
                
                for j=1:numel(actives)
                    TargetSites = connmat{actives(j)}; %#ok<*PFBNS>
                    %roll the die for each outgoing connection:
                    temp = rand(1,connects(actives(j)));
                    %collect the respective weights:
                    x = weightmat(actives(j),TargetSites);
                    temp = (x>temp).*x;
                    temp = TargetSites(temp~=0); %these are the descendants
                    
                    delaySD = delSDs(actives(j),temp);
                    at = t+sig_del(actives(j),temp);  %activation time
                    
                    for k=1:numel(temp)
                        at(k) = at(k)+randi(2*delaySD(k)+1)-(delaySD(k)+1);
                    end
                    at(at<=t) = t+1;
                    
                    temp = temp';
                    at = at';
                    
                    DrivenSites = vertcat(DrivenSites,[temp,at]);
                    pPairs = vertcat(pPairs,[preActs(j)*ones(size(at)),temp,at]);
                    [~,id] = lastwarn;
                    warning('off',id)
                end
                
                %move the clock:
                t = t+1;
                %remove past drives
                DrivenSites = DrivenSites(DrivenSites(:,2)>=t,:);
                
                if isempty(DrivenSites)
                    cShapes{n,1} = activations;
                    cWebs{n,1} = cWeb;
                    break
                else
                    %update refractory periods
                    prevZ(prevZ>1) = prevZ(prevZ>1)+1;
                    prevZ(prevZ>Taurs) = 0;
                    prevZ(actives) = 2;
                    
                    %Find driven nodes at time t and remove refractory nodes:
                    actives = DrivenSites(DrivenSites(:,2)==t);
                    actives = setdiff(actives,find(prevZ>1));
                    
                    while isempty(actives) && ~isempty(DrivenSites)
                        tD = min(DrivenSites(:,2)); %time of next driven activation
                        dt = tD-t;
                        t = tD;
                        actives = DrivenSites(DrivenSites(:,2)==t);
                        
                        %Increment node states appropriately
                        prevZ(prevZ>0) = prevZ(prevZ>0)+dt;
                        prevZ(prevZ>Taurs) = 0;
                        actives = setdiff(actives,find(prevZ>1));
                    end
                    
                    if isempty(actives)
                        cShapes{n,1} = activations;
                        cWebs{n,1} = cWeb;
                        break
                    else
                        prevZ(actives) = 1;
                        temp = [actives,t*ones(size(actives))];
                        activations = vertcat(activations,temp);
                        [~,postActs] = ismember(temp,activations,'rows');
                        
                        [~,pPair] = ismember(temp,pPairs(:,2:3),'rows');
                        cWeb = vertcat(cWeb,[pPairs(pPair,1),postActs]);
                        
                        index = true(1,size(pPairs,1));
                        index(pPair) = false;
                        pPairs = pPairs(index,:);
                        
                        cSteps = cSteps+1;
                        if cSteps>MaxSteps
                            cShapes{n,1} = activations;
                            cWebs{n,1} = cWeb;
                            disp('Maximum c-web steps exceeded!')
                            break
                        end
                    end
                end
            end
            if mod(n/Ncw*100,10)==0
                disp(n/Ncw*100)
                %             cWebs
                %             parsave(fname,{'cShapes','cWebs'},{activations,cWeb})
            end
        end
        delete(gcp)
    elseif Ncw==1
        n = 1;
        prevZ = zeros(N,1);%zeros(N,1,pTaur);
        t = 1;
        dt = 0; %#ok<NASGU>
        
        DrivenSites = [randi(N),t];   %[driven site,activation time]
        actives = DrivenSites(DrivenSites(:,2)==t);
        prevZ(actives) = 1;
        activations = DrivenSites;
        cSteps = 1;
        cWeb = zeros(0,2);
        pPairs =  zeros(0,3);%potential causal pairs
        %1st col gives presynaptic activation
        %2nd col gives postsynaptic site
        %3rd col gives postsynaptic activation time
        
        while ~isempty(DrivenSites) && cSteps<=MaxSteps
            %acts = [91 1;92 2;93 3; 94 3; 95 3];
            
            %determine activation indeces:
            [~,preActs] = ismember([actives t*ones(size(actives))],activations,'rows');
            
            for j=1:numel(actives)
                TargetSites = connmat{actives(j)}; %#ok<*PFBNS>
                %roll the die for each outgoing connection:
                temp = rand(1,connects(actives(j)));
                %collect the respective weights:
                x = weightmat(actives(j),TargetSites);
                temp = (x>temp).*x;
                temp = TargetSites(temp~=0); %these are the descendants
                
                delaySD = delSDs(actives(j),temp);
                at = t+sig_del(actives(j),temp);  %activation time
                
                for k=1:numel(temp)
                    at(k) = at(k)+randi(2*delaySD(k)+1)-(delaySD(k)+1);
                end
                at(at<=t) = t+1;
                
                temp = temp';
                at = at';
                
                DrivenSites = vertcat(DrivenSites,[temp,at]);
                pPairs = vertcat(pPairs,[preActs(j)*ones(size(at)),temp,at]);
                [~,id] = lastwarn;
                warning('off',id)
            end
            
            %move the clock:
            t = t+1;
            %remove past drives
            DrivenSites = DrivenSites(DrivenSites(:,2)>=t,:);
            
            if isempty(DrivenSites)
                cShapes{n,1} = activations;
                cWebs{n,1} = cWeb;
                break
            else
                %update refractory periods
                prevZ(prevZ>1) = prevZ(prevZ>1)+1;
                prevZ(prevZ>Taurs) = 0;
                prevZ(actives) = 2;
                
                %Find driven nodes at time t and remove refractory nodes:
                actives = DrivenSites(DrivenSites(:,2)==t);
                actives = setdiff(actives,find(prevZ>1));
                
                while isempty(actives) && ~isempty(DrivenSites)
                    tD = min(DrivenSites(:,2)); %time of next driven activation
                    dt = tD-t;
                    t = tD;
                    actives = DrivenSites(DrivenSites(:,2)==t);
                    
                    %Increment node states appropriately
                    prevZ(prevZ>0) = prevZ(prevZ>0)+dt;
                    prevZ(prevZ>Taurs) = 0;
                    actives = setdiff(actives,find(prevZ>1));
                end
                
                if isempty(actives)
                    cShapes{n,1} = activations;
                    cWebs{n,1} = cWeb;
                    break
                else
                    prevZ(actives) = 1;
                    temp = [actives,t*ones(size(actives))];
                    activations = vertcat(activations,temp);
                    [~,postActs] = ismember(temp,activations,'rows');
                    
                    [~,pPair] = ismember(temp,pPairs(:,2:3),'rows');
                    cWeb = vertcat(cWeb,[pPairs(pPair,1),postActs]);
                    
                    index = true(1,size(pPairs,1));
                    index(pPair) = false;
                    pPairs = pPairs(index,:);
                    
                    cSteps = cSteps+1;
                    if cSteps>MaxSteps
                        cShapes{n,1} = activations;
                        cWebs{n,1} = cWeb;
                        warning('Maximum c-web steps exceeded!')
                        break
                    end
                end
            end
        end
    end
end


        
