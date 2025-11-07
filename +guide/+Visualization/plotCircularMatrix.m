function plotCircularMatrix(C, labels)
% PLOTCIRCULARMATRIX  Chord-style circular plot of Yeo-7 connectivity
% with guaranteed left–right mirror ordering when both hemispheres exist.
%
% INPUTS
%   C       NxN connectivity (can be asymmetric; magnitudes are used)
%   labels  1xN cell array, e.g. {'7Networks_1 L','7Networks_1 R',...}
%
% Notes
% - Left order is forced:  VIS L, SMN L, DAN L, VAN L, LIM L, FPN L, DMN L
% - Right order is forced: DMN R, FPN R, LIM R, VAN R, DAN R, SMN R, VIS R
% - If any networks are missing on one hemisphere, the plot still works,
%   but the mirror is only perfect for the networks present on both sides.

    % ---- basic checks
    if ~ismatrix(C) || size(C,1) ~= size(C,2)
        error('C must be a square matrix.');
    end
    if numel(labels) ~= size(C,1)
        error('labels must have one name per row/column in C.');
    end

    % ---- Yeo-7 acronyms and regex to parse labels
    acr = {'VIS','SMN','DAN','VAN','LIM','FPN','DMN'};
    % more tolerant pattern: '7Networks_5 L', '7Networks-5 L', '7Networks 5 L'
    rex_num = '7Networks[\s_\-]*(\d+)\s*([LR])$';
    rex_acr = '([A-Za-z]{2,3})\s*([LR])$';   % e.g. 'FPN L'

    N = numel(labels);
    netIdx = zeros(1,N);   % 1..7, 0 = unknown
    hemi   = repmat('?',1,N);
    clean  = labels(:)';

    % ---- map labels -> (netIdx, hemi) and to acronyms
    for i = 1:N
        s = strtrim(labels{i});
        s = regexprep(s,'\s+',' '); % normalize spaces
        t = regexp(s, rex_num, 'tokens','once');
        if ~isempty(t)
            k = str2double(t{1});
            h = t{2};
            if k>=1 && k<=7, netIdx(i)=k; end
            hemi(i) = h;
            clean{i} = sprintf('%s %s', acr{netIdx(i)}, h);
        else
            t = regexp(s, rex_acr, 'tokens','once');
            if ~isempty(t)
                k = find(strcmpi(acr, upper(t{1})),1);
                if ~isempty(k), netIdx(i)=k; end
                hemi(i) = upper(t{2});
                clean{i} = sprintf('%s %s', upper(t{1}), upper(t{2}));
            else
                % leave as-is (will still plot, but won’t join mirror logic)
                clean{i} = s;
            end
        end
    end
    labels = clean;

    % ---- build explicit mirror order
    leftOrder  = [];  rightOrder = [];
    for k = 1:7
        iL = find(netIdx==k & hemi=='L', 1, 'first');
        iR = find(netIdx==k & hemi=='R', 1, 'first');
        if ~isempty(iL), leftOrder(end+1)  = iL; end %#ok<AGROW>
        if ~isempty(iR), rightOrder(end+1) = iR; end %#ok<AGROW>
    end
    % mirror: right side reversed (7..1)
    rightOrder = fliplr(rightOrder);

    % if none parsed, just place all around and continue
    if isempty(leftOrder) && isempty(rightOrder)
        warning('Could not parse Yeo-7 labels. Plotting in original order.');
        order = 1:N;
    else
        order = [leftOrder, rightOrder];
        % add any leftover nodes (unknown labels) at the end
        leftovers = setdiff(1:N, order, 'stable');
        order = [order, leftovers];
        if ~isempty(leftOrder) && ~isempty(rightOrder) && numel(leftOrder) ~= numel(rightOrder)
            warning('Unequal L/R counts (%d vs %d). Mirror guaranteed only for networks present on both sides.', ...
                    numel(leftOrder), numel(rightOrder));
        end
    end

    % ---- reorder data
    C = C(order, order);
    labels = labels(order);
    netIdx = netIdx(order);

    % ---- normalize strengths for width
    W = abs(C);
    m = max(W(:));
    if m > 0, W = W / m; end

    % ---- angles (left half pi..2pi, right half 0..pi, leftovers follow right)
    nL = numel(leftOrder);
    nR = numel(rightOrder) + numel(setdiff(order, [leftOrder rightOrder], 'stable'));
    thetaL = linspace(pi,2*pi,max(nL,1)+1); thetaL(end) = [];
    thetaR = linspace(0,pi,max(nR,1)+1);    thetaR(end) = [];
    theta  = [thetaL(1:nL), thetaR(1:nR)];   % if some side is empty, lengths handle

    % ---- assign colors by network (fixed palette; same for L/R)
    yeo7Colors = [ ...
        0.26 0.52 0.96;   % VIS (blue-ish)
        1.00 0.41 0.16;   % SMN (orange)
        0.98 0.80 0.20;   % DAN (yellow)
        0.57 0.27 0.85;   % VAN (purple)
        0.35 0.80 0.35;   % LIM (green)
        0.20 0.70 0.90;   % FPN (cyan)
        0.85 0.20 0.20];  % DMN (red)
    cmap = zeros(N,3);
    for i = 1:N
        if netIdx(i)>=1 && netIdx(i)<=7
            cmap(i,:) = yeo7Colors(netIdx(i),:);
        else
            cmap(i,:) = [0.5 0.5 0.5]; % gray for unknowns
        end
    end

    % ---- plot
    r = 1.0;  figure('Color','w'); hold on; axis equal off;

    % nodes (ticks + labels)
    for i = 1:numel(theta)
        ang = theta(i);
        plot([0.95 1.05]*r*cos(ang), [0.95 1.05]*r*sin(ang), ...
             'Color',cmap(i,:), 'LineWidth',4);
        text(1.30*r*cos(ang), 1.30*r*sin(ang), labels{i}, ...
             'HorizontalAlignment','center', 'Clipping','off', ...
             'FontWeight','bold', 'FontSize',10, 'Color',cmap(i,:));
    end

    % ribbons
    thr = 0.05;           % visibility threshold (on normalized weights)
    base = 0.02;          % thickness scale
    for i = 1:N
        for j = (i+1):N
            w = W(i,j);
            if w <= thr, continue; end
            t1 = theta(i); t2 = theta(j);
            p1 = [r*cos(t1), r*sin(t1)];
            p2 = [r*cos(t2), r*sin(t2)];
            % inward control point for a smooth Bézier arc
            pc = 0.5*(p1+p2);
            nrm = norm(pc); if nrm < eps, nrm = 1; end
            pc = pc / nrm * (0.45*r);

            tt = linspace(0,1,60).';
            curve = (1-tt).^2.*p1 + 2*(1-tt).*tt.*pc + tt.^2.*p2;

            % lateral offset to create a ribbon with thickness ~ weight
            v = curve - pc;  n = sqrt(sum(v.^2,2));  n(n<eps) = 1;
            v = v ./ n;
            off = base*w;
            curve2 = curve + off*[-v(:,2), v(:,1)]; % rotate 90° for width

            fill([curve(:,1); flipud(curve2(:,1))], ...
                 [curve(:,2); flipud(curve2(:,2))], ...
                 0.5*(cmap(i,:)+cmap(j,:)), ...
                 'EdgeColor','none', 'FaceAlpha',0.45);
        end
    end

    title('Circular Connectivity (Yeo-7, Mirrored Hemispheres)','FontSize',12);
end
