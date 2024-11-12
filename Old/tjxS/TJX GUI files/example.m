clear all;
% A=1;
% B=[2 1 5];
% %[Lia,Locb] = ismember(A,B);
% 
% if ismember(A,B)
%     disp('Yes dude');
%     [Lia,Locb]=ismember(A,B);
%     disp(Locb);
% end

% B=[1 3 40 5 6 0 9 11];
% if ismember(40,B) 
%     [Lia,Locb]=ismember(40,B);
%     for x=1:length(Lia)
%         %id = app.out.idds(x);
%         disp(x);
%     end
% end

% % Unique values
% [~,idxu,idxc] = unique(B);
% % count unique values (use histc in <=R2014b)
% [count, ~, idxcount] = histcounts(idxc,numel(idxu));
% % Where is greater than one occurence
% idxkeep = count(idxcount)>1;
% k=find(idxkeep);
% % Extract from C
% A=B(idxkeep);

% if ismember(3,B)
%     [Lia,Locb]=ismember(3,B);
%     [~,idxu,idxb] = unique(B);
%     [count, ~, idxcount] = histcounts(idxb,numel(idxu));
%     idxkeep = count(idxcount)>1;
%     k=find(idxkeep);
%     if isempty(k)
%         disp('Ninguno se repite');
%     else
%         disp(['Se repite el ' num2str(B(Locb)) ' en el lugar ' num2str(k)]);
%     end
% end
idds = [2 3 49 50 20 90];
A = [8 5 10 9 22 13];
chmap=[27 21 32 3 25 19 30 5 23 24 17 7 20 26 29 9 22 28 31 11 16 1 13 15 18 14 12 10 8 6 4 2 64 62 60 58 ...
                56 54 52 48 49 51 53 55 57 59 61 63 38 40 42 45 43 33 35 47 36 34 50 44 46 41 39 37];
            
if ismember(A, chmap)
    [Lia,Locb]=ismember(A, chmap);
    [~,idxu,idxb] = unique(A);
    [count, ~, idxcount] = histcounts(idxb,numel(idxu));
    idxkeep = count(idxcount)>1;
    k=find(idxkeep);

    if isempty(k)
        app.TextArea.Value = 'One plot';
        for i=1:length(Locb)
            for j=1:length(A)
                id = [idds(j)];
                disp(id);
                %x=chmap(Locb);
                %plot([app.UIAxes8_ x],time_step,app.out.waveforms.mean{id},'Linewidth',2,'Color','k');
                %disp(['Index ' num2str(x)]);
            %end
        end
%     else
%         app.TextArea.Value = 'Two plots';
%         for x=1:length(k)
%             id = app.out.idds(k);
%             plot(app.UIAxes8_1,time_step,app.out.waveforms.mean{id},'Linewidth',2,'Color','k');
%             hold on;
%         end
     end
end






