function plot_Fig1_NEU_separate_R2009a(dataDir, outputDir)
if nargin < 1 || isempty(dataDir), dataDir = pwd; end
if nargin < 2 || isempty(outputDir), outputDir = pwd; end
if exist(outputDir,'dir') ~= 7, mkdir(outputDir); end
files = { ...
    'NEU_application_20seed_confusion_normalized_label20_GNMFLD.csv', 'Fig1a_NEU_confusion_GNMFLD_20', 'GNMFLD'; ...
    'NEU_application_20seed_confusion_normalized_label20_GOCNMF.csv', 'Fig1b_NEU_confusion_GOCNMF_20', 'GOCNMF'; ...
    'NEU_application_20seed_confusion_normalized_label20_CGC-GOCNMF.csv', 'Fig1c_NEU_confusion_CGCGOCNMF_20', 'CGC-GOCNMF'};
labels = {'Cr','In','Pa','PS','RS','Sc'};
maximumScale = 0.75;
for k=1:size(files,1)
    M = csvread(fullfile(dataDir, files{k,1}));
    maximumScale = max(maximumScale, max(M(:)));
end
for k=1:size(files,1)
    M = csvread(fullfile(dataDir, files{k,1}));
    h = figure('Color','w','Units','pixels','Position',[100 100 610 510],'Renderer','painters');
    imagesc(M,[0 maximumScale]); axis image; axis tight; box on;
    set(gca,'XTick',1:6,'YTick',1:6,'XTickLabel',labels,'YTickLabel',labels,...
        'TickDir','in','FontName','Helvetica','FontSize',9,'LineWidth',0.8);
    xlabel('Predicted category','FontSize',10); ylabel('True category','FontSize',10);
    title(files{k,3},'FontSize',11);
    hc = colorbar; ylabel(hc,'Row-normalized proportion','FontSize',9);
    threshold = 0.5*maximumScale;
    for i=1:6
        for j=1:6
            if M(i,j) > threshold, textColor = [1 1 1]; else textColor = [0 0 0]; end
            text(j,i,sprintf('%.1f',100*M(i,j)), 'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontName','Helvetica','FontSize',8,'Color',textColor);
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
fprintf('FIG1_NEU_SEPARATE_COMPLETE=1\n');
end
