

% Analysing results

SNRl1 = zeros(nb_seed,1);
SNRl21 = zeros(nb_seed,1);

for i = 1:nb_seed % noise realizations
SNRl1(i,1) = mean(cellfun(@mean,SNR_l11{1,i}));
SNRl21(i,1) = mean(cellfun(@mean,SNR_l21{1,i})) ;

end

SNR_L1 = mean(SNRl1);
SNR_L21 = mean(SNRl21);

disp([''])
disp(['---------------------------'])
disp(['---------------------------'])
disp(['Mean SNR for single channel reconstruction: ',num2str(SNR_L1),'dB'])
disp(['Mean SNR for hyperspectral cube reconstruction: ',num2str(SNR_L21),'dB'])