%Number of simulations
n_runs = 1000;
%Amplitude of multiplicative noise
noise_scale = [0 0.5];
%Max velocity (mouse val ~18)
vel_max = [1:30];
%Min velocity (mouse val ~5)
vel_min = [1:30];
%OU time constant (mouse val ? ~0.5)
taus = [0.1:0.1:1 2];




successes = nan( length(noise_scale) , length(vel_max) , length(vel_min)  , length(taus)  ); % Store fraction successful
tts = nan( length(noise_scale) , length(vel_max) , length(vel_min)  , length(taus)  ); % Store mean time to source of successful trials
pathlength = nan( length(noise_scale) , length(vel_max) , length(vel_min)  , length(taus)  ); % Store mean pathlength of successful trials



for i = 1:length(noise_scale)
    
    for j = 1:length(vel_max)
        tic
        for k = 1:length(vel_min)
            
            parfor l = 1:length(taus)
            
                if vel_max(j) >= vel_min(k)

                    successlocal = nan(n_runs,1);
                    ttslocal = nan(n_runs,1);
                    pathlocal = nan(n_runs,1);
                    %                 tic
                    for iter = 1:n_runs
                        [ t_final, ~, ~,~, ~,~,success,~,pl] = concTrackFun_smooth_comp_ou_p( 100, 0, 2,  vel_max(j)  , vel_min(k) , noise_scale(i) , 300  , taus(l) );  % vel    vel_min     noiselevel     maxt
                        successlocal(iter) = success;
                        ttslocal(iter) = t_final;
                        pathlocal(iter) = pl;

                    end

                    %
                    successes(i,j,k,l) = nanmean(successlocal);
                    tts(i,j,k,l) = nanmean(ttslocal);
                    pathlength(i,j,k,l) = nanmean(pathlocal);
                    
                end
                
                
            end
            
            
        end
        toc
        disp( [ 'Noise: ' num2str(noise_scale(i)) '. Velmax: ' num2str(vel_max(j)) '.'] )

        
    end
    save('parameter_screen_in_prog_nobin1000_taus.mat', 'i','j','tts','successes', 'pathlength')
    
end
save('parameter_screen_nobin1000_taus.mat','tts','successes', 'pathlength')

% figure
% 
% subplot(3,3,1)
% imagesc(permute(successes(1,:,:),[2 3 1]))
% title('Percent successful, no noise')
% xlabel('V_{min} cm/s')
% ylabel('V_{max} cm/s')
% colorbar
% axis square
% 
% subplot(3,3,2)
% imagesc(permute(tts(1,:,:),[2 3 1]))
% title('Mean successful time to source, no noise')
% xlabel('V_{min} cm/s')
% ylabel('V_{max} cm/s')
% colorbar
% axis square
% 
% subplot(3,3,3)
% imagesc(permute(pathlength(1,:,:),[2 3 1]))
% title('Mean successful pathlength, no noise')
% xlabel('V_{min} cm/s')
% ylabel('V_{max} cm/s')
% colorbar
% axis square
% 
% subplot(3,3,4)
% imagesc(permute(successes(6,:,:),[2 3 1]))
% title('Percent successful, noise = 0.5')
% xlabel('V_{min} cm/s')
% ylabel('V_{max} cm/s')
% colorbar
% axis square
% 
% subplot(3,3,5)
% imagesc(permute(tts(6,:,:),[2 3 1]))
% title('Mean successful time to source, noise = 0.5')
% xlabel('V_{min} cm/s')
% ylabel('V_{max} cm/s')
% colorbar
% 
% axis square
% subplot(3,3,6)
% imagesc(permute(pathlength(6,:,:),[2 3 1]))
% title('Mean successful pathlength, noise = 0.5')
% xlabel('V_{min} cm/s')
% ylabel('V_{max} cm/s')
% colorbar
% axis square
% 
% subplot(3,3,7)
% imagesc(permute(successes(end,:,:),[2 3 1]))
% title('Percent successful, noise = 1')
% xlabel('V_{min} cm/s')
% ylabel('V_{max} cm/s')
% colorbar
% axis square
% 
% subplot(3,3,8)
% imagesc(permute(tts(end,:,:),[2 3 1]))
% title('Mean successful time to source, noise = 1')
% xlabel('V_{min} cm/s')
% ylabel('V_{max} cm/s')
% colorbar
% 
% axis square
% subplot(3,3,9)
% imagesc(permute(pathlength(end,:,:),[2 3 1]))
% title('Mean successful pathlength, noise = 1')
% xlabel('V_{min} cm/s')
% ylabel('V_{max} cm/s')
% colorbar
% axis square
% 
% figure
% 
% i = 1
% subplot(11,3,1)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Fraction Successful, Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,2)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Time to Source (sec), Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,3)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Nose Pathlength (cm), Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 2
% subplot(11,3,4)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,5)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,6)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 3
% subplot(11,3,7)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,8)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,9)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 4
% subplot(11,3,10)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,11)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,12)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 5
% subplot(11,3,13)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,14)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,15)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 6
% subplot(11,3,16)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,17)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,18)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 7
% subplot(11,3,19)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,20)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,21)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 8
% subplot(11,3,22)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,23)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,24)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 9
% subplot(11,3,25)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,26)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,27)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 10
% subplot(11,3,28)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,29)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,30)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% 
% i = 11
% subplot(11,3,31)
% imagesc(permute(successes(i,:,:),[2 3 1]))
% colorbar
% caxis([0 1])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,32)
% imagesc(permute(tts(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(tts))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square
% subplot(11,3,33)
% imagesc(permute(pathlength(i,:,:),[2 3 1]))
% colorbar
% caxis([0 max(max(max(pathlength))) ])
% title(['Noise Level: ' num2str(noise_scale(i))  ])
% xlabel('V_{min}') 
% ylabel('V_{max}')
% axis square