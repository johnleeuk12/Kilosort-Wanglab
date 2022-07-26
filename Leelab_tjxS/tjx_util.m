classdef tjx_util < handle
%{
 
%}
    
    properties
    end
    
    methods
        %%
        function app = index_stim2(self,app) 
            % function only used for 2D stim
            % input variable should contain both information from the GUI
            % (2D stim) and the xbz file containing stimulus information
            
            app.stim_info.idx = [];
            
            switch lower(app.Stim2D.Value)
                case 'trill1'
                    app.stim_info.idx = [10,22];
                case 'trill2'
                    app.stim_info.idx = [10,24];
                case 'twitter'
                    app.stim_info.idx = [9,12];
                case 'phee'
                case 'fra'
                    app.stim_info.idx = [8,3];
                case 'clicks'
                    app.stim_info.idx = [8,11];
                case ''
                    
                    
                    
            end
            app.stim_info.x_list = unique(app.out.x.stimulus_ch1(:,app.stim_info.idx(1))); % usually CF
            app.stim_info.y_list = unique(app.out.x.stimulus_ch1(:,app.stim_info.idx(2)));
            
            app.Stim1Para.Items = string([0;app.stim_info.x_list]);
            app.Stim1Para.ItemsData = [0;app.stim_info.x_list];
            app.Stim2Para.Items = string([0; app.stim_info.y_list]);
            app.Stim2Para.ItemsData = [0; app.stim_info.y_list];
            
            
        end
        
        function app = index_stim(self,app)
            app.stim_info.idx = [];
            % function to be used with 1D stim
            % input variable needs to contain the xbz file containing
            % stimulus information
            
            switch app.out.x.analysis_code
                case {1,10,62} % tones
                    app.stim_info.idx = 8;
                    app.stim_info = self.stim_label(app.stim_info, app.out.x,'tones');
                case 100    %user Calltype list
                case 102    %User Twitter1
                case 104    %User Phee
                case 2320   %VT Phee CF
                    app.stim_info.idx = 10;
                    app.stim_info = self.stim_label(app.stim_info, app.out.x,'numeric');
                case 2367   %VT Twitter CF
                    app.stim_info.idx = 8;
                    app.stim_info = self.stim_label(app.stim_info, app.out.x,'numeric');
                case 2370   %VT Twitter IPI
                    app.stim_info.idx = 8;
                    app.stim_info = self.stim_label(app.stim_info, app.out.x,'numeric');
                case 2335   %VT Trill CF
                    app.stim_info.idx = 8;
                    app.stim_info = self.stim_label(app.stim_info, app.out.x,'numeric');
                case 2338   %VT Trill FM mod
                    app.stim_info.idx = 8;
                    app.stim_info = self.stim_label(app.stim_info, app.out.x,'numeric');
                case 2341   %VT Trill AM mod
                    app.stim_info.idx = 8;
                    app.stim_info = self.stim_label(app.stim_info, app.out.x,'numeric');
                case 2343   %VT Trill AM depth
                    app.stim_info.idx = 8;
                    app.stim_info = self.stim_label(app.stim_info, app.out.x,'numeric');
                case 2339   %VT Trill FM rate
                    app.stim_info.idx = 8;
                    app.stim_info = self.stim_label(app.stim_info, app.out.x,'numeric'); 
            end
        end
        
        function stim_info = stim_label(self,stim_info,x,plot_type)
            % plot type determines the output of variables
            nreps = x.stimulus_ch1(1,4);
            TotalReps = max(x.stimulus_ch1(:,1))*nreps;
            
            switch plot_type
                case 'tones'
                    %for tuning
                    stim_info.tuning.x_label = x.stimulus_tags_ch1(stim_info.idx);
                    stim_info.tuning.x_data = x.stimulus_ch1(:,stim_info.idx);
                    stim_info.tuning.x_ticks = round(stim_info.tuning.x_data(1:4:end)*1e2)*1e-2;
                    stim_info.tuning.x_ticklabels = string(stim_info.tuning.x_ticks);
                    
                    %for raster
                    stim_info.raster.y_label = x.stimulus_tags_ch1(stim_info.idx);
                    stim_info.raster.y_ticks = [1:nreps:TotalReps]+floor(nreps);
                    stim_info.raster.y_ticklabels = round(stim_info.tuning.x_data*1e2)*1e-2;
                    
                    
                case 'numeric'
                    stim_info.tuning.x_label = x.stimulus_tags_ch1(stim_info.idx);
                    stim_info.tuning.x_data = x.stimulus_ch1(:,stim_info.idx);
                    stim_info.tuning.x_ticks = stim_info.tuning.x_data;
                    stim_info.tuning.x_ticklabels = string(stim_info.tuning.x_ticks);
                    
                    stim_info.raster.y_label = x.stimulus_tags_ch1(stim_info.idx);
                    stim_info.raster.y_ticks = [1:nreps:TotalReps]+floor(nreps);
                    stim_info.raster.y_ticklabels = stim_info.tuning.x_data;
                    
                case 'string'
                    
            end
                   
        
        
        end
    end
    
    
    
end