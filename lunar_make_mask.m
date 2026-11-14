function [mask, info] = lunar_make_mask(lat, lon, demFile, varargin)
%LUNAR_MAKE_MASK  Terrain horizon profile at a lunar south-pole site from LOLA.
%
%   mask = lunar_make_mask(lat, lon)                    % opens a file browser
%   mask = lunar_make_mask(lat, lon, 'pick')            % same
%   mask = lunar_make_mask(lat, lon, 'LDEM_80S_80MPP_ADJ.TIF')
%   [mask, info] = lunar_make_mask(lat, lon, demFile, 'Name','Value', ...)
%
%   demFile   LOLA polar elevation GeoTIFF (south polar stereographic,
%             R = 1737.4 km, true scale at the pole). Use the plain elevation
%             file, e.g. LDEM_80S_80MPP_ADJ.TIF, NOT the _HILL, _ERR, _EFFRES,
%             LDSM or LDEC versions. If the name/path is not found, the code
%             searches the current folder, this script's folder, Downloads,
%             Desktop and Documents, and finally opens a file browser.
%
%   mask      1x360, terrain elevation angle [deg] for azimuth 0..359
%             (clockwise from local north).
%
%   Options:
%     'MaxDistKm'   how far to look for terrain               (100)
%     'ObsHeightM'  observer height above ground              (0)
%     'PixelM'      pixel size [m]; default parsed from file name ("80MPP")
%     'Plot'        show the profile                          (true)
%     'SaveTo'      .mat file to store mask, lat, lon         ('')
%
%   Ray geometry was checked on a synthetic peak (azimuth matched to 0.1 deg).
%   Not run on the real LOLA files here: inspect the plotted horizon.

    if nargin < 3, demFile = 'pick'; end
    p = inputParser;
    addParameter(p,'MaxDistKm',100);
    addParameter(p,'ObsHeightM',0);
    addParameter(p,'PixelM',[]);
    addParameter(p,'Plot',true);
    addParameter(p,'SaveTo','');
    parse(p,varargin{:});
    o = p.Results;

    demFile = resolveDEM(char(demFile));
    checkName(demFile);
    fprintf('Using DEM: %s\n', demFile);

    R = 1737400;                                   % lunar reference sphere [m]
    [W,H] = rasterSize(demFile);

    if isempty(o.PixelM)
        [~,fn] = fileparts(demFile);
        tok = regexp(fn,'(\d+)MPP','tokens','once','ignorecase');
        if isempty(tok)
            error('Cannot infer pixel size from the file name. Pass ''PixelM'', e.g. 80.');
        end
        sx = str2double(tok{1});
    else
        sx = o.PixelM;
    end

    % Site in the south polar stereographic grid (x toward 90E, y toward 0E),
    % grid centered on the pole.
    lam = lon*pi/180;
    rho = 2*R*tan(pi/4 + (lat*pi/180)/2);
    xs = rho*sin(lam);  ys = rho*cos(lam);
    colS = (W+1)/2 + xs/sx;
    rowS = (H+1)/2 - ys/sx;
    if colS < 1 || colS > W || rowS < 1 || rowS > H
        error(['Site (%.3f, %.3f) lies outside this DEM (%d x %d px, %g m/px). ' ...
               'The 80S products cover only poleward of 80 S.'], lat, lon, W, H, sx);
    end

    n  = ceil(o.MaxDistKm*1000/sx) + 2;
    c1 = max(1, floor(colS)-n);  c2 = min(W, ceil(colS)+n);
    r1 = max(1, floor(rowS)-n);  r2 = min(H, ceil(rowS)+n);
    Zw = readWindow(demFile, r1, r2, c1, c2);

    zmin = min(Zw(:));  zmax = max(Zw(:));
    if zmin > 1.5e6                                % radii in meters
        Zw = Zw - R;
    elseif zmin > 1000 && zmax < 2500              % radii in km
        Zw = (Zw - 1737.4)*1000;
    end
    Zw(abs(Zw) > 2e4) = NaN;                       % nodata guard

    cs = colS - c1 + 1;  rs = rowS - r1 + 1;
    z0 = double(interp2(Zw, cs, rs, 'linear')) + o.ObsHeightM;
    if isnan(z0), error('No valid DEM height at the site.'); end

    nRay = 1440;                                   % 0.25 deg steps
    azr  = (0:nRay-1)*(360/nRay);
    a    = azr*pi/180;
    dx   = sin(lam + a);
    dy   = cos(lam + a);
    D    = sx : sx : o.MaxDistKm*1000;
    cx   = cs + dx(:)*D/sx;
    cy   = rs - dy(:)*D/sx;
    zz   = double(interp2(Zw, cx, cy, 'linear'));
    th   = D/R;
    ang  = atan2((R+zz).*cos(th) - (R+z0), (R+zz).*sin(th));
    angMax = max(ang,[],2,'omitnan')*180/pi;
    angMax(isnan(angMax)) = 0;

    bin  = mod(round(azr),360) + 1;
    mask = accumarray(bin(:), angMax, [360 1], @max).';

    info = struct('lat',lat,'lon',lon,'pixelM',sx,'siteHeightM',z0-o.ObsHeightM, ...
                  'maxMaskDeg',max(mask),'azOfMax',find(mask==max(mask),1)-1, ...
                  'window',[r1 r2 c1 c2],'file',demFile);
    fprintf('Terrain mask: site height %.0f m, max horizon %.2f deg at az %d deg.\n', ...
            info.siteHeightM, info.maxMaskDeg, info.azOfMax);

    if ~isempty(o.SaveTo), save(o.SaveTo,'mask','lat','lon'); end

    if o.Plot
        figure('Color','w','Position',[120 120 760 300]);
        area(0:359, mask, 'FaceColor',[0.55 0.55 0.55], 'EdgeColor','none'); hold on
        yline(0,'--');  xlim([0 360]);  box on
        xlabel('Azimuth (deg, clockwise from north)');  ylabel('Horizon elevation (deg)');
        title(sprintf('Terrain horizon at lat %.3f, lon %.3f (LOLA)', lat, lon));
    end
end

%% ---------------------------------------------------------------- helpers

function f = resolveDEM(f)
% Find the DEM: as given, then common folders, then a file browser.
    if ~isempty(f) && ~strcmpi(f,'pick')
        if isfile(f), return, end
        [~,nm,ex] = fileparts(f);
        names = unique({[nm ex], [nm '.TIF'], [nm '.tif']});
        home = getenv('USERPROFILE'); if isempty(home), home = getenv('HOME'); end
        folders = {pwd, fileparts(mfilename('fullpath')), ...
                   fullfile(home,'Downloads'), fullfile(home,'Desktop'), ...
                   fullfile(home,'Documents')};
        for i = 1:numel(folders)
            for j = 1:numel(names)
                c = fullfile(folders{i}, names{j});
                if isfile(c), f = c; return, end
            end
        end
        for j = 1:numel(names)                      % one recursive pass
            for i = [1 3]
                d = dir(fullfile(folders{i},'**',names{j}));
                if ~isempty(d), f = fullfile(d(1).folder, d(1).name); return, end
            end
        end
        fprintf('DEM "%s" was not found.\n  Current folder: %s\n', f, pwd);
        listLDEM(folders);
    end
    if usejava('desktop')
        [fn, pth] = uigetfile({'*.tif;*.TIF;*.tiff','GeoTIFF (*.tif)'; '*.*','All files'}, ...
                              'Select the LOLA elevation file (LDEM_80S_..._ADJ.TIF)');
        if isequal(fn,0), error('No DEM file selected.'); end
        f = fullfile(pth, fn);
    else
        error('DEM not found. Pass the full path, e.g. ''C:\\Users\\you\\Downloads\\LDEM_80S_80MPP_ADJ.TIF''.');
    end
end

function listLDEM(folders)
    found = {};
    for i = [1 3]
        d = dir(fullfile(folders{i},'**','*LDEM*'));
        for k = 1:numel(d)
            found{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
        end
    end
    if isempty(found)
        fprintf('  No *LDEM* files found in the current or Downloads folder.\n');
    else
        fprintf('  LDEM-like files found:\n');
        fprintf('    %s\n', found{:});
    end
end

function checkName(f)
    [~,fn] = fileparts(f);
    if ~isempty(regexp(fn,'(_HILL|_ERR|EFFRES|^LDSM|^LDEC)','once','ignorecase'))
        error(['"%s" is not the elevation file (it is a hillshade, error, slope or ' ...
               'count map).\nUse e.g. LDEM_80S_80MPP_ADJ.TIF.'], fn);
    end
end

function [W,H] = rasterSize(f)
    e1 = []; 
    try
        i0 = imfinfo(f);  W = i0(1).Width;  H = i0(1).Height;  return
    catch e1
    end
    try
        gi = georasterinfo(f);  H = gi.RasterSize(1);  W = gi.RasterSize(2);  return
    catch
    end
    d = dir(f);
    error(['MATLAB cannot open "%s" (%.1f MB).\n  imfinfo said: %s\n' ...
           'The download may be incomplete, or this is not a GeoTIFF. ' ...
           'LDEM_80S_80MPP_ADJ.TIF should be about 180 MB.'], f, d.bytes/1e6, e1.message);
end

function Zw = readWindow(f, r1, r2, c1, c2)
    e1 = [];  e2 = [];
    try
        raw = imread(f,'PixelRegion',{[r1 r2],[c1 c2]});
        if isinteger(raw)
            warning('DEM is integer-valued; assuming 0.5 m per count (PDS LDEM scaling).');
            Zw = single(raw)*0.5;
        else
            Zw = single(raw);
        end
        return
    catch e1
    end
    try                                         % Mapping Toolbox fallback (reads whole file)
        Z = readgeoraster(f,'OutputType','single');
        Zw = Z(r1:r2, c1:c2);
        return
    catch e2
    end
    msg2 = 'readgeoraster unavailable';  if ~isempty(e2), msg2 = e2.message; end
    error('Could not read pixels from "%s".\n  imread: %s\n  readgeoraster: %s', ...
          f, e1.message, msg2);
end
