% Inputs
leak_pam = 0.3; % Leak parameter
tn = 0.01; 
epn = 0.5;
gm_1 = 0.03; % Gamma_d1
gm_2 = 0.01; % Gamma_d2
Len = 30; % Length of system

lenG = length(gm_1);
lenG_leak = length(leak_pam);
P_ss_all = zeros(Len,lenG_leak,lenG);


%block_seq is what the determines the block size and A/B \gamma parameter,
%example 0 0 1 1 0 0 corresponds to a block system of 2 sites alternating
%ABA

% block_seq = [0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
%             0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
%             0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
%             0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
%             0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
%             0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
%             0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
%             0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
%             0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1];

block_seq = [0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 ...
    0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 ...
    0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1];


for ww=1:1:lenG
    for hh = 1:1:lenG_leak
        % NOTE: If you want to loop over gamma, change gm_2(ww) to
        % gm_1(ww), for now the functionality of the looping is limited and
        % kept this way for simplicity sake
        gm1 = gm_1(ww);
        gm2 = gm_2(ww);
        
        % calculates over different leak parameters
        gm_leak = leak_pam(hh);
        
        % This over here creates an array of gamma_d using block_seq
        for ii=1:1:Len
            if block_seq(ii) == 1
                gamma_block(ii) = gm1;
            elseif block_seq(ii) == 0
                gamma_block(ii) = gm2;
            end
        end
        
        
        % This makes an array consisting gm1+gm1 gm1+gm2 gm1+gm3 gm2+gm1...
        % for the diagonal elements of the MM matrix
        
        proper_blc = [];
        for gg=1:1:length(gamma_block)
            for rr=1:1:length(gamma_block)
                % Need to change gammma wth this somehow
                proper_blc(end+1) = -(gamma_block(gg)+gamma_block(rr));
            end
        end
        
        for Ln = 1:1:Len
            
            % Making the M matrix
            MM = zeros(Ln*Ln,Ln*Ln);
            
            for nn=1:1:Ln*Ln
                
                MM(nn,nn) = proper_blc(nn)/2;
                
                %for off-diagonal elements 
                if nn <= Ln^2-Ln
                    MM(nn,nn+Ln) = -1i*tn;
                    MM(nn+Ln,nn) = -1i*tn;
                end
                if nn <= Ln^2-1 && 0 ~= mod(nn,Ln)
                    MM(nn,nn+1) = 1i*tn;
                    MM(nn+1,nn) = 1i*tn;
                end
                
                %gm_leak correction for N contributing density matrices
                if 0 == mod(nn,Ln)
                    MM(nn,nn) = MM(nn,nn)-gm_leak/2;
                end
            end
            
            % self energies on diagonal MM zero (i.e p_11 should =0)
            for ii=1:Ln+1:Ln*Ln
                MM(ii,ii) = 0;
            end
            
            % making first row 1 0 0 0.... and set leak parameter as P_nn
            MM(1,:) = zeros(Ln*Ln,1);
            MM(1,1) = 1;
            MM(Ln^2,Ln^2) = -gm_leak;
            
            
            % Calculating the inverse matrix
            MMinv = inv(MM);
            
            % Making the v vector
            vec = zeros(Ln^2,1);
            vec(1,1) = 1;
            
            
            P_SS = MMinv*vec;
            P_ss_all(Ln,hh,ww) = P_SS(end,end)*gm_leak;
            
        end
       
    end
end

save test_data2323

clearvars -except P_ss_all MM leak_pam tn epn gm_1 gm_2 Len