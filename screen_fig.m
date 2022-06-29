%Number of simulations
n_runs = 1000;
%Amplitude of multiplicative noise
noise_scale = [0:.1:1];
%Max velocity
vel_max = [1:30];
%Min velocity
vel_min = [1:30];




successes = nan( length(noise_scale) , length(vel_max) , length(vel_min)   ); % Store fraction successful
tts = nan( length(noise_scale) , length(vel_max) , length(vel_min)   ); % Store mean time to source of successful trials
pathlength = nan( length(noise_scale) , length(vel_max) , length(vel_min)   ); % Store mean pathlength of successful trials



for i = 1:length(noise_scale)
    
    for j = 1:length(vel_max)
        tic
        parfor k = 1:length(vel_min)
            
            if vel_max(j) >= vel_min(k)
                
                successlocal = nan(n_runs,1);
                ttslocal = nan(n_runs,1);
                pathlocal = nan(n_runs,1);
                %                 tic
                for l = 1:n_runs
                    [ t_final, ~, ~,~, ~,~,success,~,pl] = concTrackFun_smooth_comp_ou_p( 100, 0, 2,  vel_max(j)  , vel_min(k) , noise_scale(i) , 300  );  % vel    vel_min     noiselevel     maxt
                    successlocal(l) = success;
                    ttslocal(l) = t_final;
                    pathlocal(l) = pl;
                    
                end
                
                %
                disp( [ 'Noise: ' num2str(noise_scale(i)) '. Velmax: ' num2str(vel_max(j)) '. Velmin: ' num2str(vel_min(k)) '.' ] )
                successes(i,j,k) = nanmean(successlocal);
                tts(i,j,k) = nanmean(ttslocal);
                pathlength(i,j,k) = nanmean(pathlocal);
                
            end
            
            
            
        end
        toc
        
    end
    save('parameter_screen_in_prog.mat', 'i','j','tts','successes', 'pathlength')
    
end
save('parameter_screen.mat','tts','successes', 'pathlength')

subplot(2,3,1)
imagesc(permute(successes(1,:,:),[2 3 1]))
title('Percent successful, no noise')
xlabel('V_{max} cm/s')
ylabel('V_{min} cm/s')
colorbar

subplot(2,3,2)
imagesc(permute(tts(1,:,:),[2 3 1]))
title('Mean successful time to source, no noise')
xlabel('V_{max} cm/s')
ylabel('V_{min} cm/s')
colorbar

subplot(2,3,3)
imagesc(permute(pathlength(1,:,:),[2 3 1]))
title('Mean successful pathlength, no noise')
xlabel('V_{max} cm/s')
ylabel('V_{min} cm/s')
colorbar


subplot(2,3,4)
imagesc(permute(successes(end,:,:),[2 3 1]))
title('Percent successful, max noise')
xlabel('V_{max} cm/s')
ylabel('V_{min} cm/s')
colorbar

subplot(2,3,5)
imagesc(permute(tts(end,:,:),[2 3 1]))
title('Mean successful time to source, max noise')
xlabel('V_{max} cm/s')
ylabel('V_{min} cm/s')
colorbar

subplot(2,3,6)
imagesc(permute(pathlength(end,:,:),[2 3 1]))
title('Mean successful pathlength, max noise')
xlabel('V_{max} cm/s')
ylabel('V_{min} cm/s')
colorbar