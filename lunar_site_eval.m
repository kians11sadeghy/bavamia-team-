function out = lunar_site_eval(lat, lon, varargin)
%LUNAR_SITE_EVAL  Sun/Earth visibility at one lunar landing site.
%
%   out = lunar_site_eval(lat, lon)
%   out = lunar_site_eval(lat, lon, 'Name','Value', ...)
%
%   lat   selenographic latitude  [deg], south negative (e.g. -89.45)
%   lon   selenographic longitude [deg], EAST positive  (e.g. 222.7 or -137.3)
%
%   Options (name-value):
%     'Name'       label for plots                         ('Landing site')
%     'Start'      UTC start, datetime or 'yyyy-MM-dd[ HH:mm]'  (today, UTC)
%     'Days'       analysis window length [days]           (30)
%     'StepMin'    time step [min]                         (60)
%     'NowOffsetH' "now" = Start + this many hours         (0)
%     'AltKm'      height above lunar reference sphere     (0)
%     'Mask'       1x360 terrain horizon [deg], az 0..359  ([] = flat 0 deg)
%     'DEM'        LOLA polar GeoTIFF path, or 'pick' for a file browser; builds the mask
%                  via lunar_make_mask.m (must be on the MATLAB path)
%     'MaxDistKm'  terrain search radius used with 'DEM'   (100)
%     'ObsHeightM' observer height above ground, with 'DEM' (0)
%     'Plot'       true/false                              (true)
%
%   Examples
%     lunar_site_eval(-86.0, 0.0);
%     lunar_site_eval(-89.45, 222.7, 'Name','Ridge', 'Start','2027-03-10', 'Days',30);
%
%   Data: JPL Horizons API (internet needed once; results are cached in
%   ./horizons_cache). Geometry: topocentric apparent azimuth/elevation of
%   the Sun (10) and Earth (399) from a surface point on the Moon. No
%   atmosphere, so the horizon is only the terrain mask you supply.
%   With the default flat mask, results are OPTIMISTIC at the south pole.
%
%   Requires MATLAB R2019b+. No toolboxes. Not validated against a
%   reference: spot-check against the Horizons web interface.

    if nargin < 2
        a = inputdlg({'Latitude (deg, south negative)','Longitude (deg east)'}, ...
                     'Landing site',1,{'-86.0','0.0'});
        if isempty(a), out = []; return, end
        lat = str2double(a{1});  lon = str2double(a{2});
        ch = questdlg('Add terrain from a LOLA elevation file (LDEM_80S_..._ADJ.TIF)?', ...
                      'Terrain','Browse for file','Flat horizon','Browse for file');
        if strcmp(ch,'Browse for file'), varargin = [varargin, {'DEM','pick'}]; end
    end
    validateattributes(lat,{'numeric'},{'scalar','real','finite','>=',-90,'<=',90});
    validateattributes(lon,{'numeric'},{'scalar','real','finite'});
    if abs(lat) > 89.999
        lat = sign(lat)*89.999;     % azimuth is undefined exactly at the pole
        warning('Latitude clamped to %.3f deg (azimuth undefined at the pole).', lat);
    end

    p = inputParser;
    addParameter(p,'Name','Landing site');
    addParameter(p,'Start',dateshift(datetime('now','TimeZone','UTC'),'start','hour'));
    addParameter(p,'Days',30);
    addParameter(p,'StepMin',60);
    addParameter(p,'NowOffsetH',0);
    addParameter(p,'AltKm',0);
    addParameter(p,'Mask',[]);
    addParameter(p,'DEM','');          % LOLA GeoTIFF -> terrain mask built for you
    addParameter(p,'MaxDistKm',100);   % terrain search radius for 'DEM'
    addParameter(p,'ObsHeightM',0);    % observer height above ground for 'DEM'
    addParameter(p,'Plot',true);
    parse(p,varargin{:});
    o = p.Results;

    cfg.startUTC         = toUTC(o.Start);
    cfg.days             = o.Days;
    cfg.stepMin          = o.StepMin;
    cfg.sunSemiDiamDeg   = 0.265;   % Sun visible if upper limb clears terrain
    cfg.earthSemiDiamDeg = 0.95;    % Earth angular semi-diameter from the Moon
    site = struct('lat',lat,'lon',lon,'alt_km',o.AltKm);
    stepH = cfg.stepMin/60;

    if isempty(o.Mask) && ~isempty(o.DEM)
        o.Mask = lunar_make_mask(lat, lon, o.DEM, 'MaxDistKm',o.MaxDistKm, ...
                                 'ObsHeightM',o.ObsHeightM, 'Plot',false);
    end
    if isempty(o.Mask)
        fprintf('Note: no terrain mask given, using a flat 0 deg horizon.\n');
    elseif numel(o.Mask) ~= 360
        error('Mask must have 360 elements (terrain elevation for az = 0..359 deg).');
    end

    %% Compute
    [t, azS, elS] = getHorizons('10',  site, cfg);
    [~, azE, elE] = getHorizons('399', site, cfg);

    sunVis   = (elS + cfg.sunSemiDiamDeg)   > interpMask(o.Mask, azS);
    earthVis = (elE + cfg.earthSemiDiamDeg) > interpMask(o.Mask, azE);

    k0 = find(t >= cfg.startUTC + hours(o.NowOffsetH), 1);
    if isempty(k0), error('"Now" falls outside the analysis window.'); end

    out.name   = o.Name;  out.lat = lat;  out.lon = lon;
    out.t = t;  out.azSun = azS;  out.elSun = elS;  out.azEarth = azE;  out.elEarth = elE;
    out.sunVisible = sunVis;  out.earthVisible = earthVis;
    out.sun   = eventMetrics(sunVis,   k0, stepH);
    out.earth = eventMetrics(earthVis, k0, stepH);
    out.sunElNow = elS(k0);
    out.sunElMax = max(elS);
    out.powerProxy = mean(max(0, sind(elS)) .* sunVis);   % flat-panel proxy

    %% Console summary
    fprintf('\n=== %s  (lat %.3f, lon %.3f) ===\n', o.Name, lat, lon);
    fprintf('Window                 : %s  +%d d, step %d min (UTC)\n', ...
            char(datetime(cfg.startUTC,'Format','yyyy-MM-dd HH:mm')), cfg.days, cfg.stepMin);
    fprintf('Sun elevation now      : %.2f deg (max in window %.2f)\n', out.sunElNow, out.sunElMax);
    fprintf('Sun                    : %s\n', eventText(out.sun));
    fprintf('Earth                  : %s\n', eventText(out.earth));
    fprintf('Sunlit fraction        : %.1f %%\n', 100*out.sun.frac);
    fprintf('Longest darkness       : %s\n', fmtDur(out.sun.longestOffH));
    fprintf('Earth visible fraction : %.1f %%\n', 100*out.earth.frac);
    fprintf('Longest Earth window   : %s\n', fmtDur(out.earth.longestOnH));
    fprintf('Flat-panel power proxy : %.3f (1 = Sun at zenith, always)\n', out.powerProxy);

    if o.Plot, makePlots(out, o, cfg, k0); end
end

%% ===================== Plotting ==========================================

function makePlots(out, o, cfg, k0)
    sunC = [0.93 0.62 0.15];  earthC = [0.22 0.54 0.87];
    N = numel(out.t);
    figure('Color','w','Position',[100 80 900 820]);
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('%s  (lat %.3f, lon %.3f)',o.Name,out.lat,out.lon),'Interpreter','none');

    % Sky view
    ax = nexttile; hold(ax,'on');
    azg = 0:359;
    mk  = interpMask(o.Mask, azg(:)).';
    fill(ax,[azg 359 0],[mk -10 -10],[0.45 0.45 0.45],'EdgeColor','none');
    w = k0 : min(k0 + round(24*60/cfg.stepMin), N);          % next 24 h
    plot(ax,out.azSun(w),  out.elSun(w),  '.','Color',sunC);
    plot(ax,out.azEarth(w),out.elEarth(w),'.','Color',earthC);
    plot(ax,out.azSun(k0),  out.elSun(k0),  'o','MarkerSize',14,'MarkerFaceColor',sunC,  'MarkerEdgeColor','k');
    plot(ax,out.azEarth(k0),out.elEarth(k0),'o','MarkerSize',10,'MarkerFaceColor',earthC,'MarkerEdgeColor','k');
    yline(ax,0,'--','Color',[0.3 0.3 0.3]);
    xlim(ax,[0 360]);  box(ax,'on');
    ylim(ax,[min(-10,min([out.elSun;out.elEarth])-2), max(20,max([out.elSun;out.elEarth])+2)]);
    xlabel(ax,'Azimuth (deg, clockwise from north)'); ylabel(ax,'Elevation (deg)');
    legend(ax,{'Terrain','Sun (24 h)','Earth (24 h)','Sun now','Earth now'},'Location','northeast');
    title(ax,'Sky view');

    % Timeline
    ax = nexttile;
    img = ones(2,N,3)*0.92;
    for c = 1:3
        r1 = img(1,:,c); r1(out.sunVisible.')   = sunC(c);   img(1,:,c) = r1;
        r2 = img(2,:,c); r2(out.earthVisible.') = earthC(c); img(2,:,c) = r2;
    end
    image(ax,[0 cfg.days],[1 2],img);
    set(ax,'YTick',[1 2],'YTickLabel',{'Sun','Earth'});
    xline(ax,hours(out.t(k0)-out.t(1))/24,'k-');
    xlim(ax,[0 cfg.days]); xlabel(ax,'Days from start');
    title(ax,'Visibility timeline (colored = above terrain horizon)');

    % Metrics
    ax = nexttile; axis(ax,'off');
    L = { sprintf('Sun elevation now : %.2f deg (max %.2f)', out.sunElNow, out.sunElMax)
          sprintf('Sun              : %s', eventText(out.sun))
          sprintf('Earth            : %s', eventText(out.earth))
          sprintf('Sunlit           : %.1f %%  (longest dark %s)', 100*out.sun.frac, fmtDur(out.sun.longestOffH))
          sprintf('Earth visible    : %.1f %%  (longest window %s)', 100*out.earth.frac, fmtDur(out.earth.longestOnH))
          sprintf('Power proxy      : %.3f', out.powerProxy) };
    text(ax,0,1,L,'VerticalAlignment','top','FontName','FixedWidth','FontSize',11);
end

%% ===================== Horizons ==========================================

function [t, az, el] = getHorizons(cmd, site, cfg)
% Topocentric apparent azimuth/elevation of body `cmd` from a lunar site.
    fmt  = 'yyyy-MM-dd HH:mm';
    stop = cfg.startUTC + days(cfg.days);
    key  = sprintf('hz2_%s_%.4f_%.4f_%.3f_%s_%d_%d', cmd, site.lat, site.lon, ...
           site.alt_km, char(datetime(cfg.startUTC,'Format','yyyyMMddHHmm')), ...
           cfg.days, cfg.stepMin);
    key   = regexprep(key,'[^\w]','_');
    cdir  = fullfile(pwd,'horizons_cache');
    cfile = fullfile(cdir,[key '.mat']);
    if isfile(cfile)
        S = load(cfile); t = S.t; az = S.az; el = S.el; return
    end

    q = @(s) ['''' s ''''];
    opts = weboptions('Timeout',90);
    r = webread('https://ssd.jpl.nasa.gov/api/horizons.api', ...
        'format','json', ...
        'COMMAND',q(cmd), 'OBJ_DATA',q('NO'), 'MAKE_EPHEM',q('YES'), ...
        'EPHEM_TYPE',q('OBSERVER'), 'CENTER',q('coord@301'), ...
        'COORD_TYPE',q('GEODETIC'), ...
        'SITE_COORD',q(sprintf('%.6f,%.6f,%.6f',mod(site.lon,360),site.lat,site.alt_km)), ...
        'START_TIME',q(char(datetime(cfg.startUTC,'Format',fmt))), ...
        'STOP_TIME', q(char(datetime(stop,'Format',fmt))), ...
        'STEP_SIZE', q(sprintf('%d min',cfg.stepMin)), ...
        'QUANTITIES',q('4'), 'ANG_FORMAT',q('DEG'), 'CSV_FORMAT',q('YES'), opts);

    txt = r.result;
    i1 = strfind(txt,'$$SOE'); i2 = strfind(txt,'$$EOE');
    if isempty(i1) || isempty(i2)
        error('Horizons returned no ephemeris. Start of reply:\n%s', txt(1:min(800,end)));
    end
    rows = splitlines(strtrim(txt(i1+5:i2-1)));
    n = numel(rows);
    t  = NaT(n,1,'TimeZone','UTC');  az = nan(n,1);  el = nan(n,1);
    for i = 1:n
        f = strsplit(rows{i},',');
        t(i) = datetime(strtrim(f{1}),'InputFormat','yyyy-MMM-dd HH:mm', ...
                        'Locale','en_US','TimeZone','UTC');
        g = strtrim(f(2:end));
        isNum = ~cellfun('isempty', regexp(g,'^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$','once'));
        v = str2double(g(isNum));                       % skips flag columns
        if numel(v) < 2
            error('Could not parse az/el from Horizons row: %s', rows{i});
        end
        az(i) = real(v(1));  el(i) = real(v(2));
    end
    if ~isfolder(cdir), mkdir(cdir); end
    save(cfile,'t','az','el');
end

%% ===================== Helpers ===========================================

function t = toUTC(s)
    if isdatetime(s), t = s; t.TimeZone = 'UTC'; return, end
    s = char(s);
    fmts = {'yyyy-MM-dd HH:mm','yyyy-MM-dd'};
    for i = 1:numel(fmts)
        try
            t = datetime(s,'InputFormat',fmts{i},'TimeZone','UTC'); return
        catch
        end
    end
    error('Start must be a datetime or text like 2027-03-10 or 2027-03-10 06:00 (UTC).');
end

function m = interpMask(mask, az)
% Periodic interpolation of a 1x360 horizon mask at azimuth(s) az [deg].
    if isempty(mask), m = zeros(size(az)); return, end
    a = [0:numel(mask)-1, numel(mask)] * (360/numel(mask));
    v = [mask(:); mask(1)].';
    m = interp1(a, v, mod(az,360));
end

function E = eventMetrics(vis, k0, stepH)
% Time to next visibility change, fraction visible, longest runs.
    vis = vis(:);  N = numel(vis);
    E.nowVisible = vis(k0);
    idx = find(vis ~= vis(k0) & (1:N).' > k0, 1);
    if isempty(idx), E.toChangeH = NaN; else, E.toChangeH = (idx-k0)*stepH; end
    E.frac = mean(vis);
    E.longestOnH  = longestRun(vis)  * stepH;
    E.longestOffH = longestRun(~vis) * stepH;
end

function L = longestRun(x)
    d = diff([0; x(:); 0]);
    s = find(d == 1);  e = find(d == -1);
    if isempty(s), L = 0; else, L = max(e - s); end
end

function s = eventText(E)
    if E.nowVisible, what = 'sets in'; else, what = 'rises in'; end
    s = sprintf('visible=%d, %s %s', E.nowVisible, what, fmtDur(E.toChangeH));
end

function s = fmtDur(h)
    if isnan(h),      s = 'beyond window';
    elseif h < 48,    s = sprintf('%.1f h', h);
    else,             s = sprintf('%.1f d', h/24);
    end
end
