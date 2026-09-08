function exportPaperFigure(fig, outputPath)
%EXPORTPAPERFIGURE Write a styled lossless PNG for paper figures.
if nargin < 2
    outputPath = fig;
    fig = gcf;
end

outputPath = char(outputPath);
[outputDir, ~, ext] = fileparts(outputPath);
if isempty(ext)
    outputPath = [outputPath, '.png'];
end
if ~isempty(outputDir) && ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

kp.plot.stylePaperFigure(fig);

if exist('export_fig', 'file') == 2
    export_fig(fig, outputPath, '-png', '-r300');
else
    exportgraphics(fig, outputPath, 'Resolution', 300);
end
end
