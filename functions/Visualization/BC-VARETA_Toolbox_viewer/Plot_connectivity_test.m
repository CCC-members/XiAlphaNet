%%
%% Testing connectivity plot
%%


addpath('data');
addpath('functions');
addpath('templates');

disp('-->> Starting process.');

% Imput params
th = 2;
load('colormap_hot_default.mat');
surf_file = 'Z:\data3_260T\data\CCLAB_DATASETS\CHBM\CHBM_ARIOSKY\update_data\BC-V_Structure\plg_segments_events\sub-CBM00001_Cl_eyes\surf\surf.mat';
connect_file = 'estimate_ec190.mat';
connect=v2m(bestVs{1});

cortex = load("Data\Atlas_Anatomical\tess_cortex_mid_high_8000V_fix.mat");
atlas = cortex.Atlas(cortex.iAtlas);
scouts = atlas.Scouts;
Vertices = cortex.Vertices;
Faces = cortex.Faces;

% Plotting surf
s = get(0, 'ScreenSize');
figure_name = 'Connectivity viewer';
connectivity_view = figure('Name',figure_name,'NumberTitle','off','Color','w','Position', [0 0 s(3) s(4)]);
%title(fig_title,'Color','k','FontSize',16);
hs = patch(...
    'Faces',            cortex.Faces, ...
    'Vertices',         cortex.Vertices,...
    'FaceVertexCData',  [], ...
    'FaceColor',        [0.6 0.6 0.6], ...
    'FaceAlpha',        0.1, ...
    'AlphaDataMapping', 'none', ...
    'EdgeColor',        'none', ...
    'EdgeAlpha',        0.1, ...
    'BackfaceLighting', 'lit', ...
    'AmbientStrength',  0.5, ...
    'DiffuseStrength',  0.5, ...
    'SpecularStrength', 0.2, ...
    'SpecularExponent', 1, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting',     'gouraud', ...
    'EdgeLighting',     'gouraud', ...
    'Tag',              'AnatSurface');
rotate3d on;
axis off;
% set(gcf,'color','w');
colormap(cmap);
colorbar;
minvalue = min(connect(:));
maxvalue = max(connect(:));
caxis([minvalue maxvalue]);

% Proting Seeds
for i=1:length(scouts)
    line(Vertices(scouts(i).Seed,1),Vertices(scouts(i).Seed,2),Vertices(scouts(i).Seed,3),...
        'LineStyle', '-', 'Marker', 'o',  'MarkerFaceColor', scouts(i).Color, 'MarkerSize', 10);
    text(Vertices(scouts(i).Seed,1),Vertices(scouts(i).Seed,2),Vertices(scouts(i).Seed,3),scouts(i).Label,'Color', scouts(i).Color);
    hold on
end

% plotting lines
for i=1:length(connect(1,:))
    for j=1:length(connect(:,1))
        w = connect(i,j);
        w = abs(w);
        if(i~=j && w>th)
            u = Vertices(scouts(i).Seed,:);
            v = Vertices(scouts(j).Seed,:);
            % Ploting Seeds
            line(u(1),u(2),u(3), 'LineStyle', '-', 'Marker', 'o',  'MarkerFaceColor', scouts(i).Color, 'MarkerSize', 10);
            line(v(1),v(2),v(3), 'LineStyle', '-', 'Marker', 'o',  'MarkerFaceColor', scouts(j).Color, 'MarkerSize', 10);
            
            %             % Define parameters of the arc.
            %             r = norm(u - v)/2;
            %             % Define the angle theta as going from 30 to 150 degrees in 100 steps.
            %             theta = linspace(0, 180, 100);
            %             % Define x and y using "Degrees" version of sin and cos.
            %             xVector = r * cosd(theta) + xCenter;
            %             yVector = r * sind(theta) + yCenter;
            %             zVector = r * sind(theta) + zCenter;
            %             % Now plot the points.
            %             LineColour=cmap(1,:);
            %             plot3(xVector,yVector,zVector,'-', 'LineWidth', w*100, 'Color', 'r');
            
            % Getting arrow distance 
            Center = [(u(1) + v(1))/2 (u(2) + v(2))/2 (u(3) + v(3))/2];
            Line = [u;Center;v];
            color = interp1(linspace(min(connect(:)),max(connect(:)),length(cmap)),cmap,connect(i,j));
            plot3(Line(:,1),Line(:,2),Line(:,3),'LineWidth', w, 'Color', color);
            num_segs = 4;
            while true
                xvals = linspace(Center(1), v(1), num_segs);
                yvals = linspace(Center(2), v(2), num_segs);
                zvals = linspace(Center(3), v(3), num_segs);
                pts = [xvals(:), yvals(:), zvals(:)];
                if(norm(Center - pts(2,:))<= 0.0030)
                    break;
                end
                num_segs = num_segs + 1;
            end
            % Plotting arrow
            if(w >= maxvalue/2)
                mArrow3(Center,pts(3,:),'stemWidth',w,'tipWidth',w,'color',color);
            else
                mArrow3(Center,pts(2,:),'stemWidth',w,'tipWidth',w,'color',color);
            end
            hold on
        end
    end
end

savefig(connectivity_view,fullfile('Outputs',[fig_title,'.fig']));
% OptionZ.FrameRate=25;
% OptionZ.Duration=20;
% OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], [pwd,filesep,'connectivity_view'],OptionZ);
% toc
% disp('done');

% v = VideoWriter('connectivity_view.avi');


disp('-->> Process finished.')
% load('colormap_hot.mat');
% 
% rng('default');
% % x = rand(64);
% data = estimate_ec190;
% thresh = 0;
% % data(data >  thresh) = 1;
% % data(data <= thresh) = 0;
% %'Colormap',myColorMap,'Label',myLabel
% cmap(191:end,:) = []; 
% circularGraph(data,'Colormap',cmap);
% colormap(cmap);
% caxis([min(data(:)) max(data(:))])
% colorbar;
