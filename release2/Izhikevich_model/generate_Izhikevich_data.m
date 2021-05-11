%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% May 2021, by Francesca Puppo
%
%%%%% Generate Izhikevich neuronal networks
%


close all
clear all

% rng(10)

% Graph parameters
Bidi = 0;
Ne_vec = [10 20 50];    %size of graph (number of neurons)
con = 2;  %based on the graph generation, the average connectivity of the network is 2*con
% for connectivity > 2 the model parameters should be changed because the activity is too high and the model isn't good to test connectivity reconstruction algorithm

% number of graphs to generate
netN = 20;
numcon_vec = zeros(1,netN);
for ineu = 3%:length(Ne_vec)
    for runid = 1:netN
        
        disp(['................Network model ' num2str(runid) '.................'])
        
        %% ---------------------------------------------------------------------
        % Parameters of simulations
        Ne=Ne_vec(ineu);
        T=200000*Ne; tau=.02; % time span and step (ms)
        fe=600; %connection strength
        fi=100;
        Ni=0;
        re=1*(1-2*rand(Ne,1));
        ri=.5+.2*(1-2*rand(Ni,1));
        a=[0.02*(ones(Ne,1)+rand(Ne,1)); 0.02+0.08*ri];
        b=[0.2*ones(Ne,1)+.01*rand(Ne,1); 0.25-0.05*ri];
        c=[-65+15*re.^2; -65*ones(Ni,1)];
        d=[6-2*re.^2; 2*ones(Ni,1)];
        
        Ntot = Ne*(Ne-1);
        
        
        %% ---------------------------------------------------------------------
        % Generate random graph
        
        aa = [repmat((1:Ne).',con,1) randi([1 Ne],[con*Ne 1])];
        aa = unique(aa,'rows');
        if Bidi==0
            aa = aa(~ismember(aa,aa(:,[2 1]),'rows'),:);                %remove bi-directional edges
        end
        aa = aa(aa(:,1)~=aa(:,2),:);                                    %remove self connecting nodes
        while ~isempty(setdiff(1:Ne,aa(:)))
            aa = [repmat((1:Ne).',con,1) randi([1 Ne],[con*Ne 1])];
            aa = unique(aa,'rows');
            if Bidi==0
                aa = aa(~ismember(aa,aa(:,[2 1]),'rows'),:);
            end
            aa = aa(aa(:,1)~=aa(:,2),:);
        end
                
        
        A = sparse(aa(:,1),aa(:,2),[.5],Ne+Ni,Ne+Ni);
        S=full(fe*A-fi*sparse([],[],[0.5],Ne+Ni,Ne+Ni));
        %      view(biograph(S.'));
        
        
        %% ---------------------------------------------------------------------
        % Simulation
        v=-70*ones(Ne+Ni,1); % Initial values of v
        u=b.*v; % Initial values of u
        firings=[]; % spike timings
        v_vec  = zeros(Ne+Ni,T+1);
        u_vec  = zeros(Ne+Ni,T+1);
        v_vec(:,1)  = v;
        u_vec(:,1)  = u;
        I=0*[10*randn(1,1);2*randn(Ni,1)]; %initialize ext current at first time step
        for t=1:10000 % simulation of 1000 ms
            if mod(t,50)==0
                I=0*[10*[randn(Ne,1).*([1;zeros(Ne-1,1)])];2*randn(Ni,1)]; % input current each 10 time steps
            end
            fired=find(v>=30); % indices of spikes
            v(fired)=c(fired);
            u(fired)=u(fired)+d(fired);
            I=I+0*sum(S(:,fired),2);
            v=v+tau*0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
            v=v+tau*0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
            u=u+tau*a.*(b.*v-u); % stability
            v_vec(:,t+1) = v;
            u_vec(:,t+1)  = u;
            
        end
                
        v_vec(:,1)  = v;
        u_vec(:,1)  = u;
        I=0;
        
        TT1 = 0;
        TT2 = round(rand*T/Ne);
        nsource_vec = randperm(Ne,Ne);
        for j=1:10
            nsource_vec = [nsource_vec randperm(Ne,Ne)];
        end
        nsource = nsource_vec(1);
        
        kk=1;
        while 1
            for t=1+TT1:TT2
                if mod(t,10)==0
                    I=[15*[randn(Ne,1).*[zeros(nsource-1,1); 1; zeros(Ne-nsource,1)]];2*randn(Ni,1)];
                    % [-5*[rand(Ne,1)];2*randn(Ni,1)]; % input current each 10 time steps
                    %     else
                    %         I=0;
                end
                fired=find(v>=30); % indices of spikes
                %             firings=[firings; t+0*fired,fired];
                v(fired)=c(fired);
                u(fired)=u(fired)+d(fired);
                I=I+sum(S(:,fired).*(1-.8*rand(Ne,length(fired))),2); %modulation of fe
                v=v+tau*0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
                v=v+tau*0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
                u=u+tau*a.*(b.*v-u); % stability
                v_vec(:,t+1) = v;
                u_vec(:,t+1)  = u;
            end
            TT1 = TT2;
            TT2 = min(round(rand*T/Ne)+TT1,T);
            kk = kk+1;
            nsource = nsource_vec(kk);
            if TT1==T
                break
            end
        end
               
        clear u_vec
        
        
        %% -----------------------------------------------------------------
        % Filtering
        v_vec(v_vec<-50) = 0;              
        
        
        %% -----------------------------------------------------------------
        % Plotting signals
        
        %             v_vec = v_vec(v_vec>-65);
        %     figure, plot(tau*(1:T+1),v_vec(1,:));
        %     figure, plot(firings(:,1),firings(:,2),'.');
        
        
        %% -----------------------------------------------------------------
        % Replace with raster signals
        v_vec_raster = cell(1,Ne);
        for i=1:Ne
            [PKS,LOCS] = findpeaks(v_vec(i,:),'MinPeakHeight',0);
            v_vec_raster{i}(LOCS,1) = 1;
        end
        time_pad = tau*(1:T+1);
        spk_wave_raster = cellfun(@(x) sparse(x),v_vec_raster,'UniformOutput',false);
        %         figure, plot(spk_wave_raster{1});
        spk_wave = v_vec;
        time_pad = [time_pad(1); time_pad(end); length(time_pad)]; 
        
    end
end


