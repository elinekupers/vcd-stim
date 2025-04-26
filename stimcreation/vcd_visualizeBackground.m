function vcd_visualizeBackground(params, bckgrnd_im, save_dir)
%% VCD function to visualize the generate background images;
%   vcd_visualizeBackground(params, bckgrnd_im, varargin)

% check inputs
if params.store_imgs && isempty(save_dir)
    fullfile(vcd_rootPath,'figs',params.disp.name,'background','visual_check');
    if ~exist(save_dir,'dir'); mkdir(save_dir); end
end

makeprettyfigures;
    
fH = figure(100);
set(fH, 'Position', [0 0 params.disp.w_pix params.disp.h_pix], 'color','w')

for jj = 1:size(bckgrnd_im,4)
    clf;
    imagesc(bckgrnd_im(:,:,:,jj)); colormap gray
    axis image
    title(sprintf('Background im %s %s %03d',gaptype, borderwidth, jj))
    set(gca, 'TickDir','out', 'LineWidth',2,'FontSize',20)
    drawnow;
    if params.store_imgs
        filename = sprintf('%03d_vcd_background_%s%s.png',jj,gaptype, borderwidth);
        imwrite(bckgrnd_im(:,:,:,jj), fullfile(vcd_rootPath,'figs',dispname,'background',filename));
        print(fH,'-dpng','-r300',fullfile(save_dir,filename));
    end
end