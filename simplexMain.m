% version April 19, 2015
% project main file - simplex
addpath('functions');
addpath('data');
clc;
%% Initialization
originalfilenames = {'lpi\_itest6.mat';'lpi\_chemcom.mat'; ...
    'lp\_agg3.mat';'lp\_80bau3b.mat'};
originalfiles = {'lpi_itest6.mat';'lpi_chemcom.mat'; ...
    'lp_agg3.mat';'lp_80bau3b.mat'};

% get files from GUI: includeFiles is the only input from the gui.
% it is a vector of 1s and 0s
totalFiles = sum(includeFiles);
filenames = cell(1,totalFiles);
files = cell(1,totalFiles);
d = 0;
for i = 1:totalFiles
    d = d+1;
    while includeFiles(d) == 0
        d = d + 1;
    end
    filenames(i) = originalfilenames(d);
    files(i) = originalfiles(d);
    
end

%%
timeRecords = zeros(numel(files),2);
problemDifficulty = zeros(numel(files),2);
d = 0;
if totalFiles ~=0
    for d = 1:totalFiles
        %% create file to store results for the current data set
        dataSetName = files{d};
        [~,tmpString,~] = fileparts(dataSetName);
        fileName = strcat(tmpString,'_simplex_results','.txt');
        createFile(fileName);
        diary on; % start writing to file
        %% import data
        fprintf('******************************************************\n \n');
        %     dataset = strcat('data/',dataSetName);
        [A,b,f] = importData(dataSetName);
        fprintf('Data set %d, %s \n', d, dataSetName);
        fprintf('Date and Time: %s \n \n', datestr(now));
        
        %% manipulate data to satisfy Ax = b and x >= 0
        [m,n] = size(A);
        lo = zeros(length(f),1);
        
        %% MATLAB linprog results on using interior-point algorithm
        fprintf('MATLAB linprog started with Interior Point\n \n');
        options = optimset('Display','iter');
        tic;
        [xM, fvalM, ~, outputM] = linprog(f,A,b,[],[],lo,[],[],options);
        timeM = toc;
        
        %% MATLAB linprog with simplex
        fprintf('MATLAB linprog started with Simplex\n \n');
        options1 = optimset('Display','iter','LargeScale','off','Simplex', 'on','TolFun',1e-8);
        tic;
        [xMS, fvalMS, ~, outputMS] = linprog(f,A,b,[],[],lo,[],[],options1);
        timeMS = toc;
        
        %% Phase 1
        fprintf('Phase 1 started \n \n');
        tic;
        [itersPH1, B_ids, fvalPH1]= phase1Setup(A,b);
        timePH1 = toc;
        if fvalPH1 >=0
            
            %% Set up phase 2
            Ahat = [A,eye(m)];
            fhat = [f;zeros(m,1)];
            
            all_ids = 1:size(Ahat,2);
            N_ids = setdiff(all_ids,B_ids);
            
            %% Phase 2
            fprintf('Phase 2 started \n \n');
            tic;
            [itersPH2, fval, x, globalSingular] = SimplexPhase2(Ahat,b,fhat,N_ids,B_ids);
            timePH2 = toc;
            
            %% Print Results
            fprintf('\nItem \t MATLAB linprog(IP) \t MATLAB Simplex \t Our Simplex\n');
            fprintf('---------------------------------------------------------------------\n');
            fprintf('fval \t %e \t \t %e \t \t %e \n', fvalM, fvalMS, fval);
            fprintf('iters \t %12d \t \t %12d \t \t %12d \n',outputM.iterations,outputMS.iterations, itersPH2);
            fprintf('time \t %e \t \t %e \t \t %e \n\n', timeM, timeMS, timePH1+timePH2);
            diary off;
            timeRecords(d,1) = timeMS;
            timeRecords(d,2) = timePH1 + timePH2;
            problemDifficulty(d,1) = itersPH2 + itersPH1;
            problemDifficulty(d,2) = globalSingular;
            %         clearvars -except files timeRecords problemDifficulty filenames d
        else
            fprintf('No feasible solution');
        end
    end
    %% Plotting Results
    close all
    bar(timeRecords);
    set(gca,'XTickLabel',filenames,'FontSize',12,'FontWeight','bold');
    if totalFiles ~= 1
        rotateXLabels(gca,90);
        legend('time taken by Matlab Simplex','time taken by our Simplex',...
            'Location','northwest');
    end
    ylabel('time in seconds','FontSize',12,'FontWeight','bold');
    title('Time comparisons','FontSize',12,'FontWeight','bold');
    colormap summer
    saveas(gca,'image1.eps','eps2c');
    %%
    figure;
    bar(problemDifficulty);
    set(gca,'XTickLabel',filenames,'FontSize',12,'FontWeight','bold');
    if totalFiles ~= 1
        rotateXLabels(gca,90);
        legend('number of iterations taken','number of singular B encountered',...
            'Location','northwest');
    end
    ylabel('count','FontSize',12,'FontWeight','bold');
    title('Analyzing the difficulty of the problems','FontSize',12,...
        'FontWeight','bold');
    saveas(gca,'image2.eps','eps2c');
    %% Create folder and move result files to that
    folderName = 'SimplexResults';
    createFolder(folderName);
    movefile('*.txt',folderName);
    movefile('*.eps',folderName);
    close all;
    input_message = 'Open the folder SimplexResults to view the results';
else
    input_message = 'Sorry! You have to select at least on dataset';
end