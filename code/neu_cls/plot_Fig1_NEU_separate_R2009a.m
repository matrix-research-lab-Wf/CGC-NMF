function plot_Fig1_NEU_separate_R2009a(dataDir, outputDir)
%PLOT_FIG1_NEU_SEPARATE_R2009A
% Generate four separate, consistently styled row-normalized confusion
% matrices for GNMFLD, ERDNMF, GOCNMF, and CGC-GOCNMF on NEU-CLS.
%
% Required CSV files in dataDir:
%   NEU_application_20seed_confusion_normalized_label20_GNMFLD.csv
%   ERDNMF_NEU_formal_20seed_confusion_normalized_label20_ERDNMF.csv
%   NEU_application_20seed_confusion_normalized_label20_GOCNMF.csv
%   NEU_application_20seed_confusion_normalized_label20_CGC-GOCNMF.csv
%
% Output:
%   Fig1a_NEU_confusion_GNMFLD_20.*
%   Fig1b_NEU_confusion_ERDNMF_20.*
%   Fig1c_NEU_confusion_GOCNMF_20.*
%   Fig1d_NEU_confusion_CGCGOCNMF_20.*
%
% MATLAB R2009a compatible.

if nargin < 1 || isempty(dataDir)
    dataDir = pwd;
end
if nargin < 2 || isempty(outputDir)
    outputDir = pwd;
end
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

files = { ...
    'NEU_application_20seed_confusion_normalized_label20_GNMFLD.csv', ...
    'Fig1a_NEU_confusion_GNMFLD_20', ...
    'GNMFLD'; ...
    'ERDNMF_NEU_formal_20seed_confusion_normalized_label20_ERDNMF.csv', ...
    'Fig1b_NEU_confusion_ERDNMF_20', ...
    'ERDNMF'; ...
    'NEU_application_20seed_confusion_normalized_label20_GOCNMF.csv', ...
    'Fig1c_NEU_confusion_GOCNMF_20', ...
    'GOCNMF'; ...
    'NEU_application_20seed_confusion_normalized_label20_CGC-GOCNMF.csv', ...
    'Fig1d_NEU_confusion_CGCGOCNMF_20', ...
    'CGC-GOCNMF'};

labels = {'Cr','In','Pa','PS','RS','Sc'};
numberOfClasses = length(labels);

% Use one common color scale for all four matrices.
maximumScale = 0.75;

for k = 1:size(files,1)
    filePath = fullfile(dataDir, files{k,1});
    if exist(filePath,'file') ~= 2
        error('Required confusion-matrix file not found: %s', filePath);
    end

    M = csvread(filePath);

    if any(size(M) ~= [numberOfClasses numberOfClasses])
        error('Confusion matrix must be %d by %d: %s', ...
            numberOfClasses, numberOfClasses, filePath);
    end
    if any(~isfinite(M(:)))
        error('Confusion matrix contains nonfinite values: %s', filePath);
    end
    if any(M(:) < -1e-12)
        error('Confusion matrix contains negative values: %s', filePath);
    end

    maximumScale = max(maximumScale, max(M(:)));
end

% Vivid colormap compatible with MATLAB R2009a.
vividMap = jet(256);

for k = 1:size(files,1)
    M = csvread(fullfile(dataDir, files{k,1}));

    h = figure( ...
        'Color','w', ...
        'Units','pixels', ...
        'Position',[100 100 610 510], ...
        'Renderer','painters');

    imagesc(M,[0 maximumScale]);
    colormap(vividMap);
    axis image;
    axis tight;
    box on;

    set(gca, ...
        'XTick',1:numberOfClasses, ...
        'YTick',1:numberOfClasses, ...
        'XTickLabel',labels, ...
        'YTickLabel',labels, ...
        'TickDir','in', ...
        'FontName','Helvetica', ...
        'FontSize',9, ...
        'LineWidth',0.8);

    xlabel('Predicted category','FontSize',10);
    ylabel('True category','FontSize',10);
   

    hc = colorbar;
    ylabel(hc,'Row-normalized proportion','FontSize',9);

    threshold = 0.5 * maximumScale;

    for i = 1:numberOfClasses
        for j = 1:numberOfClasses
            if M(i,j) > threshold
                textColor = [1 1 1];
            else
                textColor = [0 0 0];
            end

            text(j,i,sprintf('%.1f',100*M(i,j)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontName','Helvetica', ...
                'FontSize',8, ...
                'FontWeight','bold', ...
                'Color',textColor);
        end
    end

    set(h,'PaperPositionMode','auto');

    basePath = fullfile(outputDir, files{k,2});
    saveas(h,[basePath '.fig']);
    print(h,'-depsc2','-painters',[basePath '.eps']);
    print(h,'-dpdf','-painters',[basePath '.pdf']);
    print(h,'-dpng','-r600',[basePath '.png']);

    close(h);
end

fprintf('FIG1_NEU_FOUR_METHODS_SEPARATE_COMPLETE=1\n');
end
