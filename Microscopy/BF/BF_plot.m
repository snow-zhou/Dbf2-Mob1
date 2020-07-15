%% plot division to BF
mkdir('Division')
f = figure;
midslice = round(nslice/2+1);
for i = 1 : nframes
    imshow(data{1}{nslice*(i-1)+midslice}, [])
    hold on
    for n = 1 : ndivisions
        if i >= Division{n}.Budding && i <= Division{n}.Cytokinesis+1
            b = bwboundaries(imclose(trk_cells{i} == Division{n}.Mother | trk_cells{i} == Division{n}.Daughter, strel('square',3)));
            plot(smooth(b{1}(:,2)),smooth(b{1}(:,1)),'g','linewidth',1);
            d = [min(b{1}(:,2))-5, max(b{1}(:,1)+5)];
            text(d(1),d(2),num2str(n),'Color','green','FontSize',12);
        end
    end
    hold off
    saveas(f, ['Division/overlay_frame ' num2str(i)], 'png')
end