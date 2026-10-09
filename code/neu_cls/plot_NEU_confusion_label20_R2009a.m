function plot_NEU_confusion_label20_R2009a(dataDir, outputDir)
%PLOT_NEU_CONFUSION_LABEL20_R2009A
% Generate the two row-normalized NEU-CLS confusion matrices used in the
% paper at the 20% label rate.
%
% Input CSV files:
%   NEU_application_20seed_confusion_normalized_label20_GOCNMF.csv
%   NEU_application_20seed_confusion_normalized_label20_CGC-GOCNMF.csv
%
% Output files:
%   Fig_NEU_CM_GOCNMF_20.fig/.eps/.pdf/.png
%   Fig_NEU_CM_CGCGOCNMF_20.fig/.eps/.pdf/.png
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

gocPath = fullfile(dataDir, ...
    'NEU_application_20seed_confusion_normalized_label20_GOCNMF.csv');
cgcPath = fullfile(dataDir, ...
    'NEU_application_20seed_confusion_normalized_label20_CGC-GOCNMF.csv');

if exist(gocPath,'file') ~= 2
    error('File not found: %s',gocPath);
end
if exist(cgcPath,'file') ~= 2
    error('File not found: %s',cgcPath);
end

GOC = csvread(gocPath);
CGC = csvread(cgcPath);

if any(size(GOC) ~= [6 6]) || any(size(CGC) ~= [6 6])
    error('Each confusion matrix must be 6 by 6.');
end
if max(abs(sum(GOC,2)-1)) > 1e-10 || ...
        max(abs(sum(CGC,2)-1)) > 1e-10
    error('The confusion matrices are not row normalized.');
end

classLabels = {'Cr','In','Pa','PS','RS','Sc'};
maximumScale = 0.75;

local_draw_confusion_R2009a( ...
    GOC,classLabels,maximumScale, ...
    fullfile(outputDir,'Fig_NEU_CM_GOCNMF_20'));

local_draw_confusion_R2009a( ...
    CGC,classLabels,maximumScale, ...
    fullfile(outputDir,'Fig_NEU_CM_CGCGOCNMF_20'));

fprintf('NEU_CONFUSION_FIGURES_COMPLETE=1\n');
end


function local_draw_confusion_R2009a(M,classLabels,maximumScale,basePath)

figureHandle = figure( ...
    'Color','w', ...
    'Units','pixels', ...
    'Position',[100 100 610 510], ...
    'Renderer','painters');

imagesc(M,[0 maximumScale]);
axis image;
axis tight;
box on;

set(gca, ...
    'XTick',1:6, ...
    'YTick',1:6, ...
    'XTickLabel',classLabels, ...
    'YTickLabel',classLabels, ...
    'TickDir','in', ...
    'FontName','Helvetica', ...
    'FontSize',9, ...
    'LineWidth',0.8);

xlabel('Predicted defect category','FontSize',10);
ylabel('True defect category','FontSize',10);

colorbarHandle = colorbar;
ylabel(colorbarHandle,'Row-normalized proportion','FontSize',9);

colorLimits = caxis;
threshold = 0.5*(colorLimits(1)+colorLimits(2));

for i = 1:6
    for j = 1:6
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
            'Color',textColor);
    end
end

set(figureHandle,'PaperPositionMode','auto');
set(figureHandle,'InvertHardcopy','off');

saveas(figureHandle,[basePath '.fig']);
print(figureHandle,'-depsc2','-painters',[basePath '.eps']);
print(figureHandle,'-dpdf','-painters',[basePath '.pdf']);
print(figureHandle,'-dpng','-r600',[basePath '.png']);

close(figureHandle);
end
