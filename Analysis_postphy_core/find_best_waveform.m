function best_ch = find_best_waveform(x)



wv_min = zeros(1,64);
for ch = 1:64   
    wv_min(ch) = min(reshape(x.s_unit.waveforms.waveFormsMean(1,ch,10:40),1,[]));
    
end
[~,best_ch] = min(wv_min);