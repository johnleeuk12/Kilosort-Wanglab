        params = eval('parameters_distance');

        addpath('D:\DATA\M12E\Experiments')
        
        for i  = 1:length(M12E_unit_list.data)
            x = eval(M12E_unit_list.data{i,4});
            M12E_unit_list.data{i,6} = x.analysis_code;
            M12E_unit_list.data{i,7} = x.analysis_type;
        end
        
        
        
        
        
        dB = 40;
        
        u_list = find([M12E_unit_list.data{:,6}] == 1);
        p = 1;
        Pool = {};
        for i = 1:length(u_list)
            y = eval(M12E_unit_list.data{u_list(i),4});
            if y.stimulus_ch1(1,3) == dB
                unit_file_name = 'M12Eu000';
                unit_file_name = [unit_file_name(1:end-size(num2str(u_list(i)),2)) num2str(u_list(i)) '.mat'];
                x = load(unit_file_name);
                if length(x.s_unit.spiketimes)>10
                    if ~isempty(x.s_unit.templates)
                        Pool{p}.best_ch = x.s_unit.templates.best_ch;
                    else
                        data = zeros(64,60);
                        for ch = 1:64
                            temp = [];
                            temp(:,:) = x.s_unit.waveforms{1}(ch,:,:);
                            [Y,Sig,X] = svd(temp,'econ');
                            %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                            k = 1:3;
                            P = Y(:,k)*Sig(k,k)*X(:,k)';
                            data(ch,:) = mean(P,2).';
                        end
                        data = abs(data);
                        [~,Pool{p}.best_ch] = max(mean(data(:,10:40),2));
                    end
                    
                    Pool{p}.waveforms(:,:) = x.s_unit.waveforms{1}(Pool{p}.best_ch,:,:);
                    Pool{p}.spiketimes = x.s_unit.spiketimes;
                    Pool{p}.xb = y;
                    p = p+1;
                end
            end
        end