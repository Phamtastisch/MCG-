% -------------------------------------------------------------------------
% VORBEREITUNG: Herzzyklen aus dem Polar-Skript speichern
%   1. MCG_polar_v02 starten.
%   2. Im Fenster "Polarplots der Herzzyklen" oben auf "MAT speichern"
%      klicken.
%   3. Es entsteht eine Datei wie MCG_xyz_20261024_125300.mat. Darin
%      stehen x_mv, y_mv (Zyklen x 360 Phasen) und phase_deg.
%      Mehr braucht der Ortsfilter nicht.
%
% AUFRUF
%   MCG_Ortsfilter_1024                       % Dateiauswahl per Dialog
%   MCG_Ortsfilter_1024('MCG_xyz_....mat')    % Datei direkt angeben
%   MCG_Ortsfilter_1024('MCG_xyz_....mat', struct('sensor_count', 1024))
% -------------------------------------------------------------------------


% Hauptfunktion: lädt phasennormierte x/y-Spannungen aus einer MAT-Datei.
% cycle_file: Dateiname; options: Struktur zum Überschreiben der Vorgaben.
% figures enthält drei Fensterhandles. results ist beim Funktionsende noch
% leer; spätere GUI-Berechnungen liegen im userdata von Fenster 2 und 3.
function [figures, results] = MCG_Ortsfilter_1024(cycle_file, options)
%MCG_ORTSFILTER_1024 360-Grad-Zyklus aus MCG_polar_cycles-MAT-Konserve.
% Lauffaehig unter MATLAB (ab R2016b, empfohlen R2020a+ wegen UTF-8-Dateien)
% und weiterhin unter GNU Octave. Keine Toolbox erforderlich.
% Feldmodell und inverse raeumliche Musteranpassung aus MCG_Kreis_v03.
% Nur x/y werden im ebenen Strommodell verwendet; gespeichertes z entfaellt.
  % ================= EINSTELLUNGEN / RESET-AUSGANGSWERTE =================
  % Physikalisches Ersatzmodell: gerader Stromfaden in der x/y-Ebene.
  % Es beschreibt keine vollständige elektrische Volumenleitung im Körper.
  cfg.heart_diam_m = 0.05; % Länge des Stromfadens / m.
  cfg.r_heart_ohm = 20; % Ersatzwiderstand zur Spannung-Strom-Umrechnung / Ohm.
  cfg.current_gain = 1; % Dimensionsloser Stromverstärkungsfaktor.
  cfg.spatial_origin.spread_mm = 1; % Streubreite des Modellursprungs je Achse / mm: 30 = +/-15 mm; 0 = aus.
  cfg.sensor_z_m = 0.05; % Abstand der Sensorfläche von der Stromebene / m.
  cfg.grid_size_mm = 55; % Seitenlaenge des Sensorquadrats / mm: von -60 bis +60 mm.
  cfg.sensor_count = 256; % Gesamtzahl der Sensoren: Quadratzahl, z.B. 256, 1024 oder 4096.
  cfg.n_segments = 1; % Stromlinie: mehr Segmente verfeinern die Linie. Notwenig ?
  cfg.fit_angle_count = 90; % Grobsuche über 180 Grad; Vorzeichen liefert Gegenrichtung.
  cfg.fit_angle_tolerance_rad = 5; % Abbruchgenauigkeit der Winkelfeinsuche.
  cfg.current_circle_points = 200; % Auflösung des gezeichneten Herzumrisses.
  cfg.field_arrow_max_cm = 2; % Größte Bx/By-Pfeillänge pro Zyklus.
  cfg.remove_sensor_dc = true; % Räumlichen Mittelwert je Phase und Komponente entfernen.
  cfg.noiselevel_pt = 5; % Sigma (Gauss) bzw. Halbbreite (gleichverteilt) / pT.
  cfg.enable_reconstruction = true; % Haken abwählen: keine Feld-Strom-Spannungs-Rückrechnung.
  cfg.long_repeat_cycle_count = 900; % Ca. 15 Minuten Signal bei etwa 60 Herzzyklen/min.
  cfg.repeat_cycle_count = 60; % Gespeicherte Zyklen für diesen Lauf zyklisch wiederholen.
  cfg.use_spatial_filter = true; % Sensorfeld auf x/y-Modellmuster projizieren (ein Segment).
  cfg.sensor_average_samples = 1; % Zeitmittel je Sensor: 1 (aus), 3 oder 9 Phasen.
  cfg.sensor_average_choices = [1, 3, 9]; % passende  GUI-Auswahl.
  cfg.noise_distribution = 'gaussian';
  cfg.pause_s = 0.001;
  cfg.animate = true;
  cfg.measure_batch_time = true; % Laufzeitanteile im Batchbetrieb ausgeben.
  cfg.fast_batch = false; % Nur 178/180/182; Anzeige alle 100 Zyklen und am Ende.
  cfg.batch_mode = false; % Fenster 1/2 leeren; Polarplots und Sensor-Mittel in Fenster 3.
  cfg.initial_cycle_selection = 1; % 1: Mittel, 2: alle, ab 3: Einzelzyklus.
  cfg.polar_limit_mv = 0.5; % Aussenring bei fester Skalierung.
  cfg.polar_auto_scale = true; % Rauschwerte vollstaendig sichtbar machen.
  cfg.show_reconstructed = true; % Duenn blau.
  cfg.show_reference = true; % Soll aus cycles, rot.
  cfg.show_mean = true; % Mittel der Rueckrechnungen, gruen
  cfg.mean_line_width = 2.5; % Breite der durchgezogenen grünen Mittelwertlinie / pt.
  cfg.heatmap_scale = 'cycle'; % 'frame', 'cycle' oder 'fixed'.
  cfg.heatmap_limit_pt = 5000; % Fuer heatmap_scale='fixed'.
  cfg.figure1_pos = [40, 40, 1400, 1000];
  cfg.figure2_pos = [100, 80, 1200, 850];
  cfg.figure3_pos = [150, 40, 1400, 1000];
  % ====================================================================
  if nargin < 2, options = struct(); end
  names = fieldnames(options);
  for j = 1:numel(names)
    if ~isfield(cfg, names{j}), error('Unbekannte Option: %s', names{j}); end
    cfg.(names{j}) = options.(names{j});
  end
  if ~isnumeric(cfg.grid_size_mm) || ~isreal(cfg.grid_size_mm) || ...
      ~isscalar(cfg.grid_size_mm) || ~isfinite(cfg.grid_size_mm) || cfg.grid_size_mm <= 0
    error('Seitenlaenge des Sensorquadrats muss eine endliche Zahl > 0 in mm sein.');
  end
  if ~isnumeric(cfg.sensor_count) || ~isreal(cfg.sensor_count) || ...
      ~isscalar(cfg.sensor_count) || ~isfinite(cfg.sensor_count) || ...
      cfg.sensor_count < 4 || sqrt(cfg.sensor_count) ~= fix(sqrt(cfg.sensor_count))
    error('Gesamtzahl der Sensoren muss eine Quadratzahl >= 4 sein (z.B. 256, 1024, 4096).');
  end
  if ~isstruct(cfg.spatial_origin) || ~isscalar(cfg.spatial_origin) || ...
      ~isfield(cfg.spatial_origin, 'spread_mm')
    error('spatial_origin muss eine Struktur mit spread_mm sein.');
  end
  if ~isnumeric(cfg.spatial_origin.spread_mm) || ~isreal(cfg.spatial_origin.spread_mm) || ...
      ~isscalar(cfg.spatial_origin.spread_mm) || ~isfinite(cfg.spatial_origin.spread_mm) || ...
      cfg.spatial_origin.spread_mm < 0 || cfg.spatial_origin.spread_mm/2000 >= cfg.sensor_z_m
    error('Modellursprung: endliche Streubreite >= 0 mm und kleiner als der doppelte Sensorabstand erforderlich.');
  end
  sensors_per_axis = sqrt(cfg.sensor_count);
  if ~isscalar(cfg.long_repeat_cycle_count) || ~isfinite(cfg.long_repeat_cycle_count) || ...
      cfg.long_repeat_cycle_count < 1 || cfg.long_repeat_cycle_count ~= fix(cfg.long_repeat_cycle_count)
    error('Anzahl langer Zykluswiederholungen muss eine positive ganze Zahl sein.');
  end
  if ~isscalar(cfg.repeat_cycle_count) || ~isfinite(cfg.repeat_cycle_count) || ...
      cfg.repeat_cycle_count < 1 || cfg.repeat_cycle_count ~= fix(cfg.repeat_cycle_count)
    error('Anzahl wiederholter Zyklen muss eine positive ganze Zahl sein.');
  end
  if ~isscalar(cfg.sensor_average_samples) || ...
      ~ismember(cfg.sensor_average_samples, [1, 3, 9]) || ...
      ~isequal(cfg.sensor_average_choices, [1, 3, 9])
    error('Zeitmittel: nur 1, 3 oder 9 Samples erlaubt.');
  end
  if ~isscalar(cfg.noiselevel_pt) || ~isfinite(cfg.noiselevel_pt) || cfg.noiselevel_pt < 0
    error('Rauschen muss eine endliche Zahl >= 0 sein.');
  end
  if ~ismember(cfg.heatmap_scale, {'frame', 'cycle', 'fixed'}) || ...
      ~isscalar(cfg.polar_limit_mv) || ~isfinite(cfg.polar_limit_mv) || cfg.polar_limit_mv <= 0 || ...
      ~isscalar(cfg.heatmap_limit_pt) || ~isfinite(cfg.heatmap_limit_pt) || cfg.heatmap_limit_pt <= 0
    error('Ungueltige Polar- oder Heatmap-Skalierung.');
  end
  if nargin < 1 || isempty(cycle_file)
    [name, directory] = uigetfile('*.mat', 'Gespeicherten cycles-Variablensatz waehlen');
    if isequal(name, 0), figures = []; results = struct(); return; end
    cycle_file = fullfile(directory, name);
  end
  data = load(cycle_file, 'x_mv', 'y_mv', 'phase_deg');
  if ~isfield(data, 'x_mv') || ~isfield(data, 'y_mv') || ~isfield(data, 'phase_deg')
    error('cycles-MAT-Datei mit x_mv, y_mv und phase_deg erforderlich.');
  end
  if ~isequal(size(data.x_mv), size(data.y_mv)) || size(data.x_mv, 2) ~= 360 || ...
      isempty(data.x_mv) || ~isequal(data.phase_deg(:), (1:360).') || ...
      any(~isfinite(data.x_mv(:))) || any(~isfinite(data.y_mv(:)))
    error('x_mv/y_mv muessen endliche Zykluszeilen mit 360 Phasen 1...360 enthalten.');
  end
  phase = data.phase_deg(:);
  % Gleichmaessiges Raster inklusive beider Raender; mm in Meter umrechnen.
  x_grid = linspace(-cfg.grid_size_mm/2000, cfg.grid_size_mm/2000, sensors_per_axis);
  y_grid = x_grid;
  [x_obs, y_obs] = meshgrid(x_grid, y_grid);
  results = struct();
  busy = false;
  stop_requested = false;
  stop_at_phase = 0;
  completed_cycles = 0;
  active_phase = 1;
  view = struct();
  % Laufübergreifende Phasenmittel: nur bei Programmstart/Reset leeren.
  % Spalten: Sensorfeld vor/nach Ortsfilter, Rauschanteil vor/nach Ortsfilter.
  cumulative_cycle_count = 0;
  cumulative_sensor_mean = [];

  last_reconstructed_x = zeros(0, 360);
  last_reconstructed_y = zeros(0, 360);
  % Phasenweise Mittelung über Zykluszeilen; Transponieren liefert
  % 360 x 2 statt 2 x 360, passend zur späteren Stromberechnung.
  cycle_mean = [mean(data.x_mv, 1).', mean(data.y_mv, 1).'];
  fig1 = figure('name', 'MCG Ortsfilter 1024: Frontalvektor, Strom und Sensorfelder', ...
      'numbertitle', 'off', 'position', cfg.figure1_pos);
  fig2 = figure('name', 'MCG Ortsfilter 1024: Rueckgerechnete Spannungen', ...
      'numbertitle', 'off', 'position', cfg.figure2_pos);
  fig3 = figure('name', 'MCG Ortsfilter: Batch-Polarplots und gemittelte Sensorfelder', 'tag', 'batch_polars', ...
      'numbertitle', 'off', 'position', cfg.figure3_pos, 'visible', 'off');
  figures = [fig1, fig2, fig3];
  set(fig1, 'closerequestfcn', @quit_application);
  set(fig2, 'closerequestfcn', @quit_application);
  set(fig3, 'closerequestfcn', @quit_application);
  button_bar = uipanel('parent', fig2, 'units', 'normalized', 'position', [0, 0.76, 1, 0.18]);
  % Eigene, gut sichtbare Zeile für die zeitliche Sensor-Mittelung.
  average_panel = uipanel('parent', fig2, 'units', 'normalized', ...
      'position', [0, 0.94, 1, 0.06]);
  uicontrol(average_panel, 'style', 'text', 'string', 'Zeitmittel je Sensor:', ...
      'units', 'normalized', 'position', [0.01, 0.15, 0.17, 0.70], ...
      'horizontalalignment', 'left');
  average_menu = uicontrol(average_panel, 'style', 'popupmenu', ...
      'string', {'1 Sample (ohne Mittelung)', '3 Samples mitteln', '9 Samples mitteln'}, ...
      'tag', 'sensor_average_samples', 'units', 'normalized', ...
      'position', [0.18, 0.15, 0.20, 0.70], ...
      'value', find(cfg.sensor_average_choices == cfg.sensor_average_samples, 1));
  spatial_check = uicontrol(average_panel, 'style', 'checkbox', ...
      'string', 'Ortsfilter: Modell', 'tag', 'use_spatial_filter', ...
      'value', cfg.use_spatial_filter, 'units', 'normalized', ...
      'position', [0.39, 0.15, 0.17, 0.70]);
  reconstruction_check = uicontrol(average_panel, 'style', 'checkbox', ...
      'string', 'Rückrechnung', 'tag', 'enable_reconstruction', ...
      'value', cfg.enable_reconstruction, 'units', 'normalized', ...
      'position', [0.56, 0.15, 0.15, 0.70]);
  batch_check = uicontrol(average_panel, 'style', 'checkbox', 'string', 'Batchbetrieb', ...
      'tag', 'batch_mode', 'value', cfg.batch_mode, 'units', 'normalized', ...
      'position', [0.71, 0.15, 0.15, 0.70]);
  fast_batch_check = uicontrol(average_panel, 'style', 'checkbox', 'string', 'Fast Batch', ...
      'tag', 'fast_batch', 'value', cfg.fast_batch, 'units', 'normalized', ...
      'position', [0.86, 0.15, 0.14, 0.70]);
  uicontrol(button_bar, 'style', 'text', 'string', 'Sensorrauschen / pT', ...
      'units', 'normalized', 'position', [0.01, 0.52, 0.17, 0.2]);
  noise_edit = uicontrol(button_bar, 'style', 'edit', 'string', num2str(cfg.noiselevel_pt), ...
      'tag', 'sensor_noise', 'units', 'normalized', 'position', [0.18, 0.52, 0.08, 0.2]);
  noise_mode = uicontrol(button_bar, 'style', 'popupmenu', ...
      'string', {'Gauss: Sigma', 'Gleichverteilt: +/- Grenze'}, ...
      'units', 'normalized', 'position', [0.28, 0.52, 0.24, 0.2]);
  if strcmp(cfg.noise_distribution, 'uniform'), set(noise_mode, 'value', 2); end
  cycle_choices = [{'Mittlerer Herzzyklus', 'Alle Herzzyklen'}, arrayfun(@(j) sprintf('Zyklus %d', j), ...
      1:size(data.x_mv, 1), 'UniformOutput', false)];
  short_repeat_selection = numel(cycle_choices)+1;
  cycle_choices{end+1} = sprintf('%d Zyklen (wiederholt)', cfg.repeat_cycle_count);
  long_repeat_selection = numel(cycle_choices)+1;
  cycle_choices{end+1} = sprintf('%d Zyklen (ca. 15 min)', cfg.long_repeat_cycle_count);
  cycle_menu = uicontrol(button_bar, 'style', 'popupmenu', 'string', cycle_choices, 'tag', 'cycle_selection', ...
      'units', 'normalized', 'position', [0.48, 0.78, 0.29, 0.18], ...
      'value', min(cfg.initial_cycle_selection, numel(cycle_choices)));
  run_button = uicontrol(button_bar, 'style', 'pushbutton', 'string', 'Start', ...
      'tag', 'run_cycle', 'units', 'normalized', 'position', [0.01, 0.78, 0.14, 0.18], ...
      'callback', @run_cycle);
  stop_button = uicontrol(button_bar, 'style', 'pushbutton', 'string', 'Stopp', ...
      'tag', 'stop_cycle', 'enable', 'off', 'units', 'normalized', ...
      'position', [0.17, 0.78, 0.14, 0.18], 'callback', @stop_run);
  go180_button = uicontrol(button_bar, 'style', 'pushbutton', 'string', 'Go 180', ...
      'tag', 'go_180', 'units', 'normalized', 'position', [0.32, 0.78, 0.14, 0.18], ...
      'callback', @go_180);
  uicontrol(button_bar, 'style', 'pushbutton', 'string', 'Quit', ...
      'tag', 'quit_application', 'units', 'normalized', ...
      'position', [0.79, 0.78, 0.19, 0.18], 'callback', @quit_application);
  blue_check = uicontrol(button_bar, 'style', 'checkbox', 'string', 'Blau: Zyklen', ...
      'value', cfg.show_reconstructed, 'tag', 'show_blue', 'units', 'normalized', ...
      'position', [0.01, 0.27, 0.17, 0.22], 'callback', @redraw_polars);
  red_check = uicontrol(button_bar, 'style', 'checkbox', 'string', 'Rot: Soll', ...
      'value', cfg.show_reference, 'units', 'normalized', ...
      'position', [0.19, 0.27, 0.16, 0.22], 'callback', @redraw_polars);
  green_check = uicontrol(button_bar, 'style', 'checkbox', 'string', 'Gruen: Mittel', ...
      'value', cfg.show_mean, 'tag', 'show_green', 'units', 'normalized', ...
      'position', [0.36, 0.27, 0.17, 0.22], 'callback', @redraw_polars);
  auto_check = uicontrol(button_bar, 'style', 'checkbox', 'string', 'Polar: Autoskalierung', ...
      'value', cfg.polar_auto_scale, 'units', 'normalized', ...
      'position', [0.54, 0.27, 0.23, 0.22], 'callback', @redraw_polars);
  uicontrol(button_bar, 'style', 'pushbutton', 'string', 'Reset Ausgangswerte', ...
      'tag', 'reset_defaults', 'units', 'normalized', 'position', [0.79, 0.26, 0.19, 0.23], ...
      'callback', @reset_defaults);
  status = uicontrol(button_bar, 'style', 'text', 'string', 'Initialisierung', ...
      'units', 'normalized', 'position', [0.54, 0.01, 0.44, 0.22], 'tag', 'run_status', ...
      'horizontalalignment', 'left');
  phase_label = uicontrol(button_bar, 'style', 'text', 'string', 'Phase 1 / 360', ...
      'units', 'normalized', 'position', [0.01, 0.01, 0.13, 0.22]);
  phase_slider = uicontrol(button_bar, 'style', 'slider', 'min', 1, 'max', 360, ...
      'value', 1, 'sliderstep', [1/359, 10/359], 'tag', 'phase_slider', ...
      'enable', 'off', 'units', 'normalized', 'position', [0.15, 0.03, 0.37, 0.17], ...
      'callback', @scrub_phase);
  % Einheitliche Schrift fuer alle Bedienelemente und Statusangaben.
  set(findall(fig2, 'type', 'uicontrol'), ...
      'fontname', get(0, 'defaultuicontrolfontname'), 'fontunits', 'points', ...
      'fontsize', 12, 'fontweight', 'normal', 'fontangle', 'normal');
  redraw_polars();
  set(status, 'string', 'Bereit. Auswahl und Rauschen einstellen, dann Start.');

  % GUI-Callback: fordert einen kontrollierten Abbruch an.
  % Die Rechenschleifen prüfen das Flag; vollständige Ergebnisse bleiben erhalten.
  function stop_run(src, evt)
    stop_requested = true;
    if ishandle(stop_button), set(stop_button, 'enable', 'off'); end
    if ishandle(status), set(status, 'string', 'Gestoppt. Aktuelle Phase bleibt sichtbar.'); end
  end

  function go_180(src, evt)
    stop_at_phase = 180;
    if busy
      if isfield(view, 'arrow') && active_phase >= 180
        display_phase(180);
        stop_run();
      else
        set(status, 'string', 'Go 180: Lauf stoppt bei Phase 180.');
      end
    elseif isfield(view, 'arrow') || ~isempty(last_reconstructed_x)
      display_phase(180);
      stop_requested = true;
      set(status, 'string', 'Gestoppt bei Phase 180 / 360.');
      drawnow;
    else
      run_cycle(go180_button, []);
    end
  end

  % GUI-Callback: setzt das Abbruchflag und schließt die zugehörigen Fenster.
  % src und evt sind die vom Grafiksystem übergebenen Callback-Argumente.
  function quit_application(src, evt)
    stop_requested = true;
    if ishandle(fig1), delete(fig1); end
    if ishandle(fig2), delete(fig2); end
    if ishandle(fig3), delete(fig3); end
  end

  % Liefert true bei Stoppanforderung oder geschlossenen Hauptfenstern.
  % Wird auch während der inversen Rechnung als Abbruchfunktion verwendet.
  function value = should_stop()
    value = stop_requested || ~ishandle(fig1) || ~ishandle(fig2);
  end

  % GUI-Callback: setzt Bedienelemente auf die beim Aufruf geltenden Vorgaben.
  % Löscht bisherige Ergebnisse und Achsen; während eines Laufs ohne Wirkung.
  function reset_defaults(src, evt)
    if busy, return; end
    set(noise_edit, 'string', num2str(cfg.noiselevel_pt));
    set(noise_mode, 'value', 1+strcmp(cfg.noise_distribution, 'uniform'));
    set(cycle_menu, 'value', min(cfg.initial_cycle_selection, numel(cycle_choices)));
    set(blue_check, 'value', cfg.show_reconstructed);
    set(red_check, 'value', cfg.show_reference);
    set(green_check, 'value', cfg.show_mean);
    set(auto_check, 'value', cfg.polar_auto_scale);
    set(batch_check, 'value', cfg.batch_mode);
    set(fast_batch_check, 'value', cfg.fast_batch);
    set(average_menu, 'value', find(cfg.sensor_average_choices == cfg.sensor_average_samples, 1));
    set(spatial_check, 'value', cfg.use_spatial_filter);
    set(reconstruction_check, 'value', cfg.enable_reconstruction);
    last_reconstructed_x = zeros(0, 360); last_reconstructed_y = zeros(0, 360);
    completed_cycles = 0;
    cumulative_cycle_count = 0;
    cumulative_sensor_mean = [];
    view = struct(); active_phase = 1; stop_at_phase = 0;
    set(phase_slider, 'value', 1, 'enable', 'off');
    set(phase_label, 'string', 'Phase 1 / 360');
    results = struct();
    set(fig2, 'userdata', results);
    delete(findall(fig1, 'type', 'axes'));
    delete(findall(fig2, 'type', 'axes'));
    delete(findall(fig3, 'type', 'axes'));
    set(fig3, 'visible', 'off');
    redraw_polars();
    set(status, 'string', 'Ausgangswerte wiederhergestellt. Bereit fuer Start.');
  end

  % GUI-Callback: zeichnet x/y-Zyklen mit den aktuellen Anzeigeoptionen.
  % Wählt das Batchfenster oder Hauptfenster und bestimmt eine gemeinsame Skala.
  % Verändert nur die Anzeige, nicht die bereits berechneten Spannungen.
  function redraw_polars(src, evt)
    if ~ishandle(fig1), return; end
    % Drei Einzelphasen ergeben keinen vollstaendigen Polarzyklus.
    if ~isempty(last_reconstructed_x) && size(last_reconstructed_x, 2) ~= 360, return; end
    if get(batch_check, 'value')
      target_fig = fig3;
      delete(findall(target_fig, 'tag', 'voltage_polar_axes'));
      if isempty(last_reconstructed_x), return; end
      set(target_fig, 'visible', 'on');
      position_x = [0.12, 0.72, 0.28, 0.24];
      position_y = [0.60, 0.72, 0.28, 0.24];
    else
      target_fig = fig1;
      delete(findall(target_fig, 'tag', 'voltage_polar_axes'));
      position_x = [0.10, 0.07, 0.30, 0.23];
      position_y = [0.59, 0.07, 0.30, 0.23];
    end
    plot_cfg = cfg;
    plot_cfg.show_reconstructed = logical(get(blue_check, 'value'));
    plot_cfg.show_reference = logical(get(red_check, 'value'));
    plot_cfg.show_mean = logical(get(green_check, 'value'));
    plot_cfg.polar_auto_scale = logical(get(auto_check, 'value'));
    plot_cfg.polar_limit_mv = cfg.polar_limit_mv;
    if plot_cfg.polar_auto_scale
      % Größter Spannungsbetrag plus 5 % Rand; mindestens der Vorgabewert.
      plot_cfg.polar_limit_mv = max(cfg.polar_limit_mv, ...
          1.05*max(abs([cycle_mean(:); last_reconstructed_x(:); last_reconstructed_y(:)])));
    end
    cycle_voltage_polar(target_fig, position_x, cycle_mean(:, 1), ...
        sprintf('Rueckgerechnetes x: %d Zyklen', size(last_reconstructed_x, 1)), last_reconstructed_x, plot_cfg);
    cycle_voltage_polar(target_fig, position_y, cycle_mean(:, 2), ...
        sprintf('Rueckgerechnetes y: %d Zyklen', size(last_reconstructed_y, 1)), last_reconstructed_y, plot_cfg);
    update_phase_points();
    drawnow;
  end

  % GUI-Callback: führt den ausgewählten mittleren oder einzelne Zyklen aus.
  % Ablauf: Spannung -> Strom -> Sensorfeld -> Rauschen -> inverse Anpassung
  % -> Spannung. Speichert vollständige Zyklen und optional Batchlaufzeiten.
  % Die Referenzdaten dienen der Simulation und Fehlerprüfung, nicht dem Fit.
  function run_cycle(src, evt)
    if busy || ~ishandle(fig1) || ~ishandle(fig2), return; end
    level = str2double(get(noise_edit, 'string'));
    if ~isfinite(level) || level < 0
      set(status, 'string', 'Rauschen: endliche Zahl >= 0 eingeben.'); return;
    end
    if get(spatial_check, 'value') && cfg.n_segments ~= 1
      set(status, 'string', 'Linearer Ortsfilter benötigt n_segments=1.');
      return;
    end
    if isequal(src, run_button), stop_at_phase = 0; end
    busy = true;
    stop_requested = false;
    fast_batch = logical(get(fast_batch_check, 'value'));
    batch_mode = fast_batch || logical(get(batch_check, 'value'));
    computed_phases = 1:360;
    if fast_batch, computed_phases = [178, 180, 182]; end
    n_phases = numel(computed_phases);
    % Wechsel des Phasenrasters darf keine unterschiedlichen Mittelwerte mischen.
    if ~isempty(cumulative_sensor_mean) && size(cumulative_sensor_mean, 1) ~= numel(x_obs)*n_phases*3
      cumulative_cycle_count = 0;
      cumulative_sensor_mean = [];
    end
    measure_time = batch_mode && cfg.measure_batch_time;
    % Wandzeit einschließlich GUI-Ereignisverarbeitung; keine reine CPU-Zeit.
    run_timer = tic;
    timing_rows = zeros(0, 4); % Vorwärts, invers, Polarplots, Gesamt / s.
    if batch_mode
      delete(findall(fig1, 'type', 'axes'));
      delete(findall(fig2, 'type', 'axes'));
      delete(findall(fig3, 'type', 'axes'));
      set(fig3, 'visible', 'on');
    else
      delete(findall(fig3, 'type', 'axes'));
      set(fig3, 'visible', 'off');
    end
    view = struct(); active_phase = 1;
    set(phase_slider, 'enable', 'off', 'value', 1);
    set(run_button, 'enable', 'off');
    set(stop_button, 'enable', 'on');
    set(batch_check, 'enable', 'off');
    set(fast_batch_check, 'enable', 'off');
    set(average_menu, 'enable', 'off');
    set(spatial_check, 'enable', 'off');
    set(reconstruction_check, 'enable', 'off');
    % Lokales onCleanup-Objekt: wird bei jedem return/Fehler zerstoert und
    % ruft finish_run auf (kein 'clear' noetig, in MATLAB-Nested-Functions heikel).
    cleanup = onCleanup(@finish_run); %#ok<NASGU>
    local_cfg = cfg;
    local_cfg.noiselevel_pt = level;
    local_cfg.enable_reconstruction = logical(get(reconstruction_check, 'value'));
    if ~local_cfg.enable_reconstruction
      delete(findall(fig2, 'type', 'axes'));
    end
    local_cfg.use_spatial_filter = logical(get(spatial_check, 'value'));
    average_selection = get(average_menu, 'value');
    local_cfg.sensor_average_samples = cfg.sensor_average_choices(average_selection);
    local_cfg.fast_batch = fast_batch;
    generated_phases = computed_phases;
    if fast_batch
      half_window = (local_cfg.sensor_average_samples-1)/2;
      % Alle benoetigten Nachbarn einmal berechnen; ueberlappende Fenster teilen Daten.
      generated_phases = unique(reshape(bsxfun(@plus, computed_phases(:), ...
          -half_window:half_window), 1, []));
    end
    [unused, output_indices] = ismember(computed_phases, generated_phases);
    local_cfg.should_stop = @should_stop;
    if get(noise_mode, 'value') == 1, local_cfg.noise_distribution = 'gaussian';
    else, local_cfg.noise_distribution = 'uniform'; end
    choice = get(cycle_menu, 'value');
    cycle_mean = [mean(data.x_mv, 1).', mean(data.y_mv, 1).'];
    if choice == long_repeat_selection
      source_cycles = mod(0:cfg.long_repeat_cycle_count-1, size(data.x_mv, 1))+1;
    elseif choice == short_repeat_selection
      % Wiederholt die vorhandenen Zyklen; Rauschen wird pro Durchlauf neu gezogen.
      source_cycles = mod(0:cfg.repeat_cycle_count-1, size(data.x_mv, 1))+1;
    elseif choice == 2, source_cycles = 1:size(data.x_mv, 1);
    elseif choice == 1, source_cycles = 0;
    else, source_cycles = choice-2; end
    sensor_rms_squared = zeros(0, 4); % Feld vorher/nachher, Rauschen vorher/nachher / T².
    reconstructed_x = zeros(0, n_phases);
    reconstructed_y = zeros(0, n_phases);
    last_reconstructed_x = reconstructed_x; last_reconstructed_y = reconstructed_y;
    completed_cycles = 0;
    for batch_cycle = 1:numel(source_cycles)
    drawnow;
    if should_stop(), return; end
    update_display = ~fast_batch || mod(batch_cycle, 100) == 0 || batch_cycle == numel(source_cycles);
    cycle_timer = tic;
    forward_timer = tic;
    source_cycle = source_cycles(batch_cycle);
    if source_cycle == 0
      p_mv = cycle_mean;
    else
      p_mv = [data.x_mv(source_cycle, :).', data.y_mv(source_cycle, :).'];
    end
    generated_p_mv = p_mv(generated_phases, :);
    p_mv = p_mv(computed_phases, :);
    % cycles-Konvention: aVF=-y. Keine Spiegelung im Feldmodell.
    original_leads = frontal_to_goldberger(p_mv);
    % Ohmsches Ersatzmodell: mV / Ohm ergibt mA; Faktor 1000 ergibt µA.
    current_ua = local_cfg.current_gain*1000/local_cfg.r_heart_ohm*generated_p_mv;
    % Ein gemeinsamer Modellursprung je Herzzyklus; die Sensorplatte bleibt fest.
    % Das wahre Feld entsteht am Ursprung 0, Filter/Fit nehmen den gezogenen Ursprung an.
    half_spread_m = local_cfg.spatial_origin.spread_mm/2000;
    local_cfg.model_origin_m = half_spread_m*(2*rand(1, 3)-1);
    simulation_cfg = local_cfg;
    simulation_cfg.model_origin_m = [0, 0, 0];
    clean = zeros(numel(y_grid), numel(x_grid), numel(generated_phases), 3);
    for k = 1:numel(generated_phases)
      drawnow;
      if should_stop(), return; end
      [bx, by, bz] = current_line_field(current_ua(k, 1), current_ua(k, 2), x_obs, y_obs, simulation_cfg);
      clean(:, :, k, 1) = bx; clean(:, :, k, 2) = by; clean(:, :, k, 3) = bz;
    end
    % 1 pT = 1e-12 T: Rauschpegel vor Addition in SI-Einheiten umrechnen.
    measured = clean + sensor_noise(size(clean), level*1e-12, local_cfg.noise_distribution);
    % Sensordaten und inverse Modelltemplates identisch raeumlich zentrieren.
    if local_cfg.remove_sensor_dc
      % Mittel über Dimension 1 und 2: räumlicher Offset je Phase und
      % Feldkomponente. bsxfun verteilt den Offset über alle Sensoren.
      measured = bsxfun(@minus, measured, mean(mean(measured, 1), 2));
      clean = bsxfun(@minus, clean, mean(mean(clean, 1), 2));
    end
    % Zeitfilter erst nach Erzeugung und Zentrierung der Sensordaten.
    % Dimension 3 ist die Zyklusphase; Sensorpositionen und B-Komponenten
    % werden niemals miteinander gemischt. Rohdaten zum Vergleich aufbewahren.
    raw_sensor_values_t = measured;
    raw_clean_sensor_values_t = clean;
    measured = temporal_sensor_average(measured, local_cfg.sensor_average_samples);
    clean = temporal_sensor_average(clean, local_cfg.sensor_average_samples);
    % Nur die drei Zielphasen weiterverarbeiten. Ihre Fenster liegen komplett
    % im berechneten Bereich; dessen zyklische Raender werden nicht verwendet.
    measured = measured(:, :, output_indices, :);
    clean = clean(:, :, output_indices, :);
    raw_sensor_values_t = raw_sensor_values_t(:, :, output_indices, :);
    raw_clean_sensor_values_t = raw_clean_sensor_values_t(:, :, output_indices, :);
    current_ua = current_ua(output_indices, :);
    % Modellfilter nach Zeitmittelung: zwei physikalische Feldmuster
    % erhalten, dazu orthogonale Störanteile aus den Sensorwerten entfernen.
    % Vor dem Filter liegende Daten für einen direkten Vergleich speichern.
    pre_spatial_sensor_values_t = measured;
    spatial_filter = struct('enabled', local_cfg.use_spatial_filter);
    if local_cfg.use_spatial_filter
      [measured, spatial_filter] = spatial_sensor_projection(measured, x_obs, y_obs, local_cfg);
    end
    drawnow;
    if should_stop(), return; end
    forward_s = toc(forward_timer); % Einschließlich Rauschen und Zentrierung.
    if ~batch_mode
    delete(findall(fig2, 'type', 'axes'));
    set(0, 'currentfigure', fig1);
    delete(findall(fig1, 'type', 'axes'));
    set(fig1, 'colormap', jet(256));
    ax_vector = axes('parent', fig1, 'position', [0.06, 0.72, 0.24, 0.22]);
    plot(ax_vector, p_mv(:, 1), p_mv(:, 2), '-', 'color', [0.75, 0.75, 0.75]); hold(ax_vector, 'on');
    arrow = quiver(ax_vector, 0, 0, p_mv(1, 1), p_mv(1, 2), 0, 'r', 'linewidth', 2);
    point = plot(ax_vector, p_mv(1, 1), p_mv(1, 2), 'ro');
    vector_limit = max(0.05, 1.1*max(abs(p_mv(:))));
    axis(ax_vector, 'equal'); xlim(ax_vector, [-1, 1]*vector_limit); ylim(ax_vector, [-1, 1]*vector_limit);
    grid(ax_vector, 'on'); xlabel(ax_vector, 'x / mV'); ylabel(ax_vector, 'y / mV (-aVF)');
    ax_current = axes('parent', fig1, 'position', [0.39, 0.72, 0.24, 0.22]);
    phi = linspace(0, 2*pi, cfg.current_circle_points); radius = local_cfg.heart_diam_m/2;
    plot(ax_current, radius*cos(phi), radius*sin(phi), 'k-'); hold(ax_current, 'on');
    stream = plot(ax_current, [0, 0], [0, 0], 'b-', 'linewidth', 2);
    tip = plot(ax_current, 0, 0, 'bo');
    % Eine Stromlinie je Phase; bereits durchlaufene Phasen bleiben sichtbar.
    % Euklidischer Strombetrag sqrt(Ix²+Iy²), numerisch mit hypot berechnet.
    amplitudes = hypot(current_ua(:, 1), current_ua(:, 2));
    colors = jet(256);
    amplitude_max = max([amplitudes; 1e-12]);
    history_lines = empty_graphics(360, 1);
    for sample = 1:360
      if mod(sample, 10) == 1
        drawnow;
        if should_stop(), return; end
      end
      if amplitudes(sample) > 0
        % Strom durch seinen Betrag teilen: Einheitsrichtung. Die Länge
        % bleibt der Herzradius; Stromstärke erscheint über die Linienfarbe.
        endpoint_history = radius*current_ua(sample, :)/amplitudes(sample);
      else
        endpoint_history = [0, 0];
      end
      % Linear von [0, Maximalstrom] auf die 256 Farbindizes [1,256] abbilden.
      color_index = 1+round(255*amplitudes(sample)/amplitude_max);
      history_lines(sample) = plot(ax_current, ...
          [-endpoint_history(1), endpoint_history(1)], ...
          [-endpoint_history(2), endpoint_history(2)], ...
          'color', colors(color_index, :), 'linewidth', 0.8, 'visible', 'off');
    end
    caxis(ax_current, [0, amplitude_max]);
    current_colorbar = colorbar(ax_current);
    if isprop(current_colorbar, 'Label')
      set(get(current_colorbar, 'Label'), 'String', '|I| / uA'); % MATLAB
    else
      ylabel(current_colorbar, '|I| / uA'); % Octave
    end
    axis(ax_current, 'equal'); xlim(ax_current, 1.15*radius*[-1, 1]); ylim(ax_current, 1.15*radius*[-1, 1]);
    grid(ax_current, 'on'); xlabel(ax_current, 'x / m'); ylabel(ax_current, 'y / m');
    map_handles = empty_graphics(1, 3); map_axes = empty_graphics(1, 3); components = {'Bx', 'By', 'Bz'};
    % Tesla -> pT; symmetrische Farbskala um null, gemeinsames Zyklusmaximum.
    field_limit = max(1, max(abs(measured(:)))*1e12);
    for component = 1:3
      drawnow;
      if should_stop(), return; end
      map_axes(component) = axes('parent', fig1, 'position', [0.055+(component-1)*0.33, 0.40, 0.25, 0.23]);
      map_handles(component) = imagesc(map_axes(component), x_grid*100, y_grid*100, measured(:, :, 1, component)*1e12);
      axis(map_axes(component), 'image'); set(map_axes(component), 'ydir', 'normal');
      caxis(map_axes(component), [-1, 1]*field_limit); colorbar(map_axes(component));
      set(map_axes(component), 'tag', ['heatmap_' components{component}]);
      title(map_axes(component), [components{component} ' Sensorwerte / pT']);
      xlabel(map_axes(component), 'x / cm'); ylabel(map_axes(component), 'y / cm');
    end
    ax_quiver = axes('parent', fig1, 'position', [0.72, 0.72, 0.24, 0.22]);
    % Feste Pfeilskalierung fuer den gesamten Zyklus, keine Normierung je Frame.
    % Betrag des Felds in der Sensorebene / pT. Skalenfaktor / cm pro pT
    % hält die Pfeillängen während des ganzen Zyklus vergleichbar.
    xy_fields = hypot(measured(:, :, :, 1), measured(:, :, :, 2))*1e12;
    quiver_scale = cfg.field_arrow_max_cm/max(1, max(xy_fields(:))); % Maximum gemäß cfg.field_arrow_max_cm.
    field_quiver = quiver(ax_quiver, x_obs*100, y_obs*100, ...
        measured(:, :, 1, 1)*1e12*quiver_scale, ...
        measured(:, :, 1, 2)*1e12*quiver_scale, 0, 'k');
    axis(ax_quiver, 'equal'); grid(ax_quiver, 'on');
    xlim(ax_quiver, [min(x_grid)*100-2.5, max(x_grid)*100+2.5]);
    ylim(ax_quiver, [min(y_grid)*100-2.5, max(y_grid)*100+2.5]);
    xlabel(ax_quiver, 'x / cm'); ylabel(ax_quiver, 'y / cm');
    title(ax_quiver, 'Bx/By Sensorvektoren (feste Pfeilskalierung)');
    redraw_polars();
    view = struct();
    view.arrow = arrow;
    view.point = point;
    view.ax_vector = ax_vector;
    view.current_ua = current_ua;
    view.p_mv = p_mv;
    view.radius = radius;
    view.stream = stream;
    view.tip = tip;
    view.history_lines = history_lines;
    view.ax_current = ax_current;
    view.field_quiver = field_quiver;
    view.measured = measured;
    view.quiver_scale = quiver_scale;
    view.map_handles = map_handles;
    view.map_axes = map_axes;
    view.field_limit = field_limit;
    display_phase(1);
    drawnow;
    if should_stop(), return; end
    end % Keine Aktualisierung von Fenster 1 im Batchbetrieb.
    if local_cfg.enable_reconstruction
    if update_display
      set(status, 'string', sprintf('Zyklus %d/%d: Spannungen aus Sensoren rueckrechnen ...', batch_cycle, numel(source_cycles)));
    end
    drawnow;
    if should_stop(), return; end
    inverse_timer = tic;
    fit = reconstruct_spatial_field(measured, x_grid, y_grid, local_cfg);
    if should_stop() || fit.aborted, return; end
    % Umkehrung des Ohmschen Ersatzmodells: U_mV = I_µA * R / (1000*Gain).
    reconstructed_mv = fit.current_ua*local_cfg.r_heart_ohm/(1000*local_cfg.current_gain);
    reconstructed_leads = frontal_to_goldberger(reconstructed_mv);
    inverse_s = toc(inverse_timer); % Einschließlich Umrechnung in Spannungen.
    reconstructed_x(end+1, :) = reconstructed_mv(:, 1).';
    reconstructed_y(end+1, :) = reconstructed_mv(:, 2).';
    last_reconstructed_x = reconstructed_x; last_reconstructed_y = reconstructed_y;
    plot_timer = tic;
    if ~fast_batch, redraw_polars(); end
    plot_s = toc(plot_timer); % drawnow ist in redraw_polars enthalten.
    else
      % Ortsfilter und RMS-Messung bleiben aktiv. Keine inverse Winkelsuche,
      % keine Spannungsumrechnung und keine rekonstruierten Polarlinien.
      fit = struct('enabled', false);
      reconstructed_mv = zeros(0, 2);
      reconstructed_leads = zeros(0, 3);
      inverse_s = 0;
      plot_s = 0;
    end
    completed_cycles = completed_cycles+1;
    % RMS über alle Sensoren, 360 Phasen und Bx/By/Bz: sqrt(mean(B²)).
    % Feld-RMS enthält Nutzsignal plus Rauschen. Separater Rausch-RMS
    % verwendet das zeitlich genauso gemittelte, rauschfreie Simulationsfeld.
    % "Vorher" liegt nach dem Zeitmittel, aber vor dem räumlichen Filter.
    noise_before = pre_spatial_sensor_values_t-clean;
    noise_after = measured-clean;
    sensor_rms_squared(end+1, :) = [mean(pre_spatial_sensor_values_t(:).^2), ...
        mean(measured(:).^2), mean(noise_before(:).^2), mean(noise_after(:).^2)];

    % Erst Sensorwerte gleicher Phase über die Zyklen mitteln, dann RMS.
    % Online-Mittel m_neu = m_alt + (Wert-m_alt)/Anzahl benötigt keinen
    % Speicher für alte Zyklen. Sensor-/Phasen-/Komponentenreihenfolge bleibt
    % durch (:)-Stapelung identisch. Auch unterschiedliche Sollzyklen werden
    % phasenweise gemittelt; der Rauschanteil bezieht sich jeweils auf ihr Modell.
    current_sensor_vectors = [pre_spatial_sensor_values_t(:), measured(:), ...
        noise_before(:), noise_after(:)];
    cumulative_cycle_count = cumulative_cycle_count+1;
    if cumulative_cycle_count == 1
      cumulative_sensor_mean = current_sensor_vectors;
    else
      cumulative_sensor_mean = cumulative_sensor_mean + ...
          (current_sensor_vectors-cumulative_sensor_mean)/cumulative_cycle_count;
    end
    cumulative_rms_pt = sqrt(mean(cumulative_sensor_mean.^2, 1))*1e12;

    results = struct('phase_deg', phase(computed_phases), 'original_xy_mv', p_mv, ...
        'reconstructed_xy_mv', reconstructed_mv, 'original_leads_mv', original_leads, ...
        'reconstructed_leads_mv', reconstructed_leads, 'sensor_values_t', measured, ...
        'clean_sensor_values_t', clean, 'fit', fit, 'cfg', local_cfg, 'cycle_file', cycle_file);
    results.model_origin_offset_m = local_cfg.model_origin_m;
    results.true_origin_m = [0, 0, 0];
    results.actual_sensor_positions_m = cat(3, x_obs, y_obs, local_cfg.sensor_z_m*ones(size(x_obs)));
    results.sensor_mean_since_reset_t = reshape(cumulative_sensor_mean(:, 2), size(measured));
    results.batch_sample_phases_deg = [178, 180, 182];
    [unused, sample_indices] = ismember([178, 180, 182], computed_phases);
    results.batch_sample_mean_t = results.sensor_mean_since_reset_t(:, :, sample_indices, :);
    % Nach Projektion enthaelt der Vergleich mit clean auch Modell-/Lagefehler.
    results.sensor_rms.columns = {'field_before_pt', 'field_after_pt', ...
        'noise_before_pt', 'noise_after_pt'};
    results.sensor_rms.per_cycle_pt = sqrt(sensor_rms_squared)*1e12;
    % Gleiche Anzahl Messwerte pro Zyklus: mittlere Quadrate zusammenfassen,
    % erst danach die Wurzel ziehen (nicht die einzelnen RMS-Werte mitteln).
    results.sensor_rms.averaged_since_reset_pt = cumulative_rms_pt;
    results.sensor_rms.cycles_since_reset = cumulative_cycle_count;
    results.sensor_rms.total_pt = sqrt(mean(sensor_rms_squared, 1))*1e12;
    results.pre_spatial_sensor_values_t = pre_spatial_sensor_values_t;
    results.spatial_filter = spatial_filter;
    results.raw_sensor_values_t = raw_sensor_values_t;
    results.raw_clean_sensor_values_t = raw_clean_sensor_values_t;
    results.sensor_average_samples = local_cfg.sensor_average_samples;
    results.generated_phases_deg = generated_phases;
    results.reconstructed_cycles_x_mv = reconstructed_x;
    results.reconstructed_cycles_y_mv = reconstructed_y;
    if isempty(reconstructed_x)
      results.reconstructed_mean_xy_mv = zeros(0, 2);
    else
      results.reconstructed_mean_xy_mv = [mean(reconstructed_x, 1).', mean(reconstructed_y, 1).'];
    end
    results.source_cycles = source_cycles(1:batch_cycle);
    set(fig2, 'userdata', results);
    if ~batch_mode && local_cfg.enable_reconstruction
    delete(findall(fig2, 'type', 'axes'));
    reference = [p_mv, original_leads]; recovered = [reconstructed_mv, reconstructed_leads];
    labels = {'x', 'y (-aVF)', 'aVR', 'aVL', 'aVF'};
    for lead = 1:5
      drawnow;
      if should_stop(), return; end
      ax = axes('parent', fig2, 'position', [0.09, 0.62-(lead-1)*0.135, 0.85, 0.105]);
      h_recovered = plot(ax, phase, recovered(:, lead), 'b-', 'linewidth', 0.8, ...
          'tag', 'recovered_voltage');
      hold(ax, 'on');
      h_reference = plot(ax, phase, reference(:, lead), 'r-', 'linewidth', 1.6, ...
          'tag', 'reference_voltage');
      hold(ax, 'off');
      grid(ax, 'on'); xlim(ax, [1, 360]); ylabel(ax, [labels{lead} ' / mV']);
      if lead == 1, title(ax, 'Spannungen nach einem Herzzyklus'); legend(ax, ...
          [h_reference, h_recovered], {'Konserve (Soll)', 'Aus verrauschten Sensoren'}, 'location', 'northeast'); end
      if lead == 5, xlabel(ax, 'Herzzyklusphase / Grad'); end
    end
    end % Batchbetrieb zeigt Polarplots und gemittelte Sensorfelder in Fenster 3.
    if ~batch_mode
      if stop_at_phase == 180 || cfg.animate, frames = 1:360; else, frames = 360; end
      for k = frames
        drawnow;
        if should_stop(), return; end
        display_phase(k);
        if stop_at_phase == 180 && k == 180
          stop_run();
          return;
        end
        set(status, 'string', sprintf('Zyklus %d/%d | Phase %d / 360 | Rauschen %.4g pT', ...
            batch_cycle, numel(source_cycles), k, level));
        drawnow;
        if should_stop(), return; end
        if cfg.animate && cfg.pause_s > 0
          pause_timer = tic;
          while toc(pause_timer) < cfg.pause_s
            pause(min(0.01, max(0, cfg.pause_s-toc(pause_timer))));
            if should_stop(), return; end
          end
        end
      end
    end
    % RMSE über beide Komponenten und alle Phasen / mV:
    % sqrt(Mittelwert((U_rekonstruiert-U_Soll)²)).
    if update_display
    if local_cfg.enable_reconstruction
      rmse = sqrt(mean((reconstructed_mv(:)-p_mv(:)).^2));
      set(status, 'string', sprintf('Zyklus %d fertig | x/y-RMSE %.5g mV', completed_cycles, rmse));
    else
      set(status, 'string', sprintf('Zyklus %d fertig | Rückrechnung aus | Sensor-RMS wird gesammelt', completed_cycles));
    end
    if batch_mode
      heatmap_timer = tic;
      redraw_batch_means();
      plot_s = plot_s + toc(heatmap_timer);
      set(fig3, 'userdata', results);
      figure(fig3);
      if stop_at_phase == 180
        display_phase(180);
        stop_run();
      end
    else
      figure(fig2);
    end
    end % Anzeige im Fast Batch nur alle 100 Zyklen und am Laufende.
    drawnow;
    if should_stop(), return; end
    if measure_time
      timing_rows(end+1, :) = [forward_s, inverse_s, plot_s, toc(cycle_timer)];
      % Nur abgeschlossene Zyklen aufnehmen. Gesamt enthält zusätzlich
      % Ergebnisverwaltung und GUI-Aufwand außerhalb der Polarzeichnung.
      results.timing.columns = {'forward_s', 'inverse_s', 'polar_plot_s', 'total_s'};
      results.timing.per_cycle_s = timing_rows;
      results.timing.sum_s = sum(timing_rows, 1);
      results.timing.run_elapsed_s = toc(run_timer);
      if update_display
      fprintf('Batch Zyklus %d: vorwärts %.3f s | Rückrechnung %.3f s | Polarplots %.3f s | gesamt %.3f s\n', ...
          completed_cycles, timing_rows(end, :));
      set(fig2, 'userdata', results);
      set(fig3, 'userdata', results);
      set(status, 'string', sprintf('Batch: %d Zyklen | vorwärts %.2f s | Rückrechnung %.2f s | Polarplots %.2f s | Laufzeit %.2f s', ...
          completed_cycles, results.timing.sum_s(1:3), results.timing.run_elapsed_s));
      end
    end
    end
    if measure_time && ~isempty(timing_rows)
      fprintf('Batch gesamt (%d Zyklen): vorwärts %.3f s | Rückrechnung %.3f s | Polarplots %.3f s | Laufzeit %.3f s\n', ...
          completed_cycles, sum(timing_rows(:, 1:3), 1), toc(run_timer));
    end
    if ~isempty(sensor_rms_squared)
      rms_pt = sqrt(mean(sensor_rms_squared, 1))*1e12;
      fprintf('Sensor-RMS gesamt (%d Zyklen, Zeitmittel %d Samples):\n', ...
          completed_cycles, local_cfg.sensor_average_samples);
      fprintf('  Sensorfeld vor/nach Ortsfilter: %.6g / %.6g pT\n', rms_pt(1:2));
      fprintf('  Abweichung vor/nach Ortsfilter: %.6g / %.6g pT (Rauschen und ggf. Modell-/Lagefehler)\n', rms_pt(3:4));
      fprintf('  RMS gemittelt seit Anfang/Reset (%d Zyklen): Sensorfeld %.6g / %.6g pT | Rauschen %.6g / %.6g pT (vor/nach Ortsfilter)\n', ...
          cumulative_cycle_count, cumulative_rms_pt);
    end
    % finish_run laeuft automatisch, wenn cleanup beim Verlassen zerstoert wird.
  end

  % Phasenweise Sensor-Mittelwerte seit Start/Reset, nach Zeit-/Ortsfilter.
  % Quiver zeigt den gemittelten Bx/By-Vektor je Sensor, keine Sollrichtung.
  function redraw_batch_means()
    samples = results.batch_sample_mean_t*1e12;
    sample_phases = results.batch_sample_phases_deg;
    components = {'Bx', 'By', 'Bz'};
    color_limit = max(1e-9, max(abs(samples(:))));
    xy_magnitude = hypot(samples(:, :, :, 1), samples(:, :, :, 2));
    arrow_scale = cfg.field_arrow_max_cm/max(1e-9, max(xy_magnitude(:)));
    stride = max(1, ceil(numel(x_grid)/12));
    ix = 1:stride:numel(x_grid); iy = 1:stride:numel(y_grid);
    set(fig3, 'colormap', jet(256));
    for sample = 1:3
      for component = 1:3
        tag = sprintf('batch_mean_%s_%d', components{component}, sample_phases(sample));
        ax = findall(fig3, 'tag', tag);
        if isempty(ax)
          ax = axes('parent', fig3, 'position', ...
              [0.065+(sample-1)*0.32, 0.49-(component-1)*0.21, 0.245, 0.145], ...
              'tag', tag);
          im = imagesc(ax, x_grid*100, y_grid*100, samples(:, :, sample, component));
          set(im, 'tag', 'batch_mean_image');
          axis(ax, 'image'); set(ax, 'ydir', 'normal');
          colorbar(ax); xlabel(ax, 'x / cm'); ylabel(ax, 'y / cm');
          if component == 3
            hold(ax, 'on');
            quiver(ax, x_obs(iy, ix)*100, y_obs(iy, ix)*100, ...
                samples(iy, ix, sample, 1)*arrow_scale, ...
                samples(iy, ix, sample, 2)*arrow_scale, 0, 'k', ...
                'linewidth', 1.2, 'tag', sprintf('batch_mean_quiver_%d', sample_phases(sample)));
            hold(ax, 'off');
          end
        else
          image_handle = findall(ax, 'type', 'image');
          set(image_handle, 'cdata', samples(:, :, sample, component));
          if component == 3
            q = findall(ax, 'tag', sprintf('batch_mean_quiver_%d', sample_phases(sample)));
            set(q, 'udata', samples(iy, ix, sample, 1)*arrow_scale, ...
                'vdata', samples(iy, ix, sample, 2)*arrow_scale);
          end
        end
        caxis(ax, [-1, 1]*color_limit);
        label = sprintf('%s / pT | Phase %d | Mittel: %d Zyklen', ...
            components{component}, sample_phases(sample), cumulative_cycle_count);
        if component == 3, label = {label; 'Pfeile: mittlere Bx/By'}; end
        title(ax, label, 'fontsize', 10);
      end
      drawnow;
      if should_stop(), return; end
    end
  end

  % Der Schieberegler verändert nur die Anzeige der vorhandenen Zyklusdaten.
  function scrub_phase(src, evt)
    if busy, return; end
    display_phase(round(get(phase_slider, 'value')));
    set(status, 'string', sprintf('Gestoppt | Phase %d / 360', active_phase));
    drawnow;
  end

  function display_phase(k)
    previous_phase = active_phase;
    active_phase = max(1, min(360, round(k)));
    k = active_phase;
    set(phase_slider, 'value', k);
    set(phase_label, 'string', sprintf('Phase %d / 360', k));
    if isfield(view, 'arrow') && ishandle(view.arrow)
      set(view.arrow, 'udata', view.p_mv(k, 1), 'vdata', view.p_mv(k, 2));
      set(view.point, 'xdata', view.p_mv(k, 1), 'ydata', view.p_mv(k, 2));
      title(view.ax_vector, sprintf('Frontalvektor | Phase %d / 360 Grad', k));
      amplitude = norm(view.current_ua(k, :));
      if amplitude > 0, endpoint = view.radius*view.current_ua(k, :)/amplitude; else, endpoint = [0, 0]; end
      set(view.stream, 'xdata', [-endpoint(1), endpoint(1)], 'ydata', [-endpoint(2), endpoint(2)]);
      set(view.tip, 'xdata', endpoint(1), 'ydata', endpoint(2));
      if k < previous_phase
        set(view.history_lines(k+1:previous_phase), 'visible', 'off');
      elseif k > previous_phase
        set(view.history_lines(previous_phase+1:k), 'visible', 'on');
      end
      set(view.history_lines(k), 'visible', 'on');
      title(view.ax_current, sprintf('%d Stromlinien | |I| = %.4g uA', k, amplitude));
      set(view.field_quiver, 'udata', view.measured(:, :, k, 1)*1e12*view.quiver_scale, ...
          'vdata', view.measured(:, :, k, 2)*1e12*view.quiver_scale);
      for component = 1:3
        set(view.map_handles(component), 'cdata', view.measured(:, :, k, component)*1e12);
        if strcmp(cfg.heatmap_scale, 'frame')
          frame_values = view.measured(:, :, k, component)*1e12;
          color_limit = max(1e-9, max(abs(frame_values(:))));
        elseif strcmp(cfg.heatmap_scale, 'fixed')
          color_limit = cfg.heatmap_limit_pt;
        else
          color_limit = view.field_limit;
        end
        caxis(view.map_axes(component), [-1, 1]*color_limit);
      end
    end
    update_phase_points();
  end

  % Punkte übernehmen ihre Koordinaten direkt aus den sichtbaren Kurven.
  function update_phase_points()
    tags = {'recovered_voltage', 'reference_voltage', 'polar_cycle', ...
        'polar_reference', 'polar_reconstructed_mean'};
    for fig = figures
      if ~ishandle(fig), continue; end
      for t = 1:numel(tags)
        curves = findall(fig, 'tag', tags{t});
        for n = 1:numel(curves)
          curve = curves(n);
          x = get(curve, 'xdata'); y = get(curve, 'ydata');
          if numel(x) < active_phase || numel(y) < active_phase, continue; end
          marker = getappdata(curve, 'phase_point');
          if isempty(marker) || ~ishandle(marker)
            marker = line('parent', get(curve, 'parent'), 'xdata', x(active_phase), ...
                'ydata', y(active_phase), 'linestyle', 'none', 'marker', 'o', ...
                'markersize', 7, 'color', get(curve, 'color'), ...
                'markerfacecolor', get(curve, 'color'), 'tag', 'phase_point', ...
                'handlevisibility', 'off');
            setappdata(curve, 'phase_point', marker);
          else
            set(marker, 'xdata', x(active_phase), 'ydata', y(active_phase));
          end
        end
      end
    end
  end

  % Aufräumfunktion für onCleanup: gibt Bedienelemente wieder frei.
  % Wird bei normalem Ende, frühem return und auch bei einem Fehler ausgeführt.
  function finish_run()
    busy = false;
    if ishandle(phase_slider) && ishandle(fig1) && ishandle(fig2)
      if isfield(view, 'arrow') || (~isempty(last_reconstructed_x) && size(last_reconstructed_x, 2) == 360)
        set(phase_slider, 'enable', 'on');
        display_phase(active_phase);
      end
    end
    if ishandle(run_button), set(run_button, 'enable', 'on'); end
    if ishandle(stop_button), set(stop_button, 'enable', 'off'); end
    if ishandle(batch_check), set(batch_check, 'enable', 'on'); end
    if ishandle(fast_batch_check), set(fast_batch_check, 'enable', 'on'); end
    if ishandle(average_menu), set(average_menu, 'enable', 'on'); end
    if ishandle(spatial_check), set(spatial_check, 'enable', 'on'); end
    if ishandle(reconstruction_check), set(reconstruction_check, 'enable', 'on'); end
    if stop_requested && ishandle(status)
      set(status, 'string', sprintf('Gestoppt bei Phase %d / 360. %d vollstaendige Zyklen berechnet. Start beginnt einen neuen Lauf.', active_phase, completed_cycles));
    end
  end
end

% Zeichnet Spannung gegen Zyklusphase als eigenen Polarplot.
% voltage: Referenzspalte / mV; cycles: Zeilen einzelner Zyklen / mV.
% cfg steuert Grenzwert und Sichtbarkeit; Rückgabe ax ist das Achsenhandle.
function ax = cycle_voltage_polar(fig, position, voltage, label, cycles, cfg)
  ax = axes('parent', fig, 'position', position, 'tag', 'voltage_polar_axes');
  hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
  xlim(ax, [-1.23, 1.23]); ylim(ax, [-1.23, 1.23]);
  circle = linspace(0, 2*pi, 361);
  patch('parent', ax, 'xdata', 0.5*cos(circle(:)), 'ydata', 0.5*sin(circle(:)), ...
      'facecolor', [0.85, 0.85, 0.85], 'edgecolor', 'none');
  limit = cfg.polar_limit_mv;
  ticks = limit*[-0.5, 0, 0.5, 1];
  for value = ticks
    r = (value+limit)/(2*limit);
    if value == 0, style = 'k--'; else, style = ':'; end
    plot(ax, r*cos(circle), r*sin(circle), style);
  end
  for angle = 0:60:300
    phi = angle*pi/180;
    plot(ax, [0, cos(phi)], [0, sin(phi)], '-', 'color', [0.7, 0.7, 0.7]);
    shown_angle = angle; if angle == 0, shown_angle = 360; end
    text(1.13*cos(phi), 1.13*sin(phi), sprintf('%d°', shown_angle), ...
        'parent', ax, 'horizontalalignment', 'center', 'fontsize', 8);
  end
  % Winkel = zeitliche Zyklusphase in rad. Spannung wird auf [0,1]
  % abgebildet: -limit -> Zentrum, 0 -> Radius 0.5, +limit -> Außenring.
  % NaN unterbricht nur die Darstellung außerhalb des Spannungsbereichs.
  theta = (1:360).'*pi/180;
  r = (voltage(:)+limit)/(2*limit);
  r(r < 0 | r > 1) = NaN;
  % Dünne blaue Einzelzyklen zuerst zeichnen.
  if cfg.show_reconstructed
    for cycle = 1:size(cycles, 1)
      r_cycle = (cycles(cycle, :).'+limit)/(2*limit);
      r_cycle(r_cycle < 0 | r_cycle > 1) = NaN;
      plot(ax, r_cycle.*cos(theta), r_cycle.*sin(theta), 'b-', 'linewidth', 0.7, ...
          'tag', 'polar_cycle');
    end
  end
  for value = ticks
    r = (value+limit)/(2*limit);
    text(r*cos(pi/4), r*sin(pi/4), sprintf('%+.2g mV', value), ...
        'parent', ax, 'horizontalalignment', 'center', 'fontsize', 8, ...
        'backgroundcolor', [1, 1, 1]);
  end
  % Rote Solllinie über den Einzelzyklen zeichnen.
  if cfg.show_reference
    r_reference = (voltage(:)+limit)/(2*limit);
    r_reference(r_reference < 0 | r_reference > 1) = NaN;
    plot(ax, r_reference.*cos(theta), r_reference.*sin(theta), 'r-', ...
        'linewidth', 2, 'tag', 'polar_reference');
  end
  % Mittelwert zuletzt zeichnen: andere Kurven dürfen die grüne Linie
  % nicht abschnittsweise überdecken. Den Übergang von Phase 360 zu 1
  % schließen, damit die Kurve am Phasenwechsel fortlaufend bleibt.
  if ~isempty(cycles) && cfg.show_mean
    mean_radius = (mean(cycles, 1).'+limit)/(2*limit);
    mean_radius(mean_radius < 0 | mean_radius > 1) = NaN;
    mean_theta = [theta; theta(1)+2*pi];
    mean_radius = [mean_radius; mean_radius(1)];
    plot(ax, mean_radius.*cos(mean_theta), mean_radius.*sin(mean_theta), ...
        'linestyle', '-', 'color', [0, 0.8, 0], ...
        'linewidth', cfg.mean_line_width, 'tag', 'polar_reconstructed_mean');
  end
  title(ax, label);
  text(0, -1.35, 'Blau: rueckgerechnet | Rot: Soll | Gruen: Mittel', ...
      'parent', ax, 'horizontalalignment', 'center', 'fontsize', 8);
end

% Projiziert Zeilenvektoren p=[x,y] auf aVR, aVL und aVF.
% Eingabe und Ausgabe haben dieselbe Spannungseinheit; drei Ausgabespalten.
% Die drei Achsenvektoren summieren sich zu null: aVR+aVL+aVF=0.
function leads = frontal_to_goldberger(p)
  % Skalarprodukte mit den Goldberger-Achsen:
  % aVR=(-sqrt(3)/2,1/2), aVL=(sqrt(3)/2,1/2), aVF=(0,-1).
  leads = [-sqrt(3)/2*p(:, 1)+0.5*p(:, 2), ...
            sqrt(3)/2*p(:, 1)+0.5*p(:, 2), -p(:, 2)];
end

% Erzeugt unabhängiges additives Rauschen in der angegebenen Arraygröße.
% level_t ist in Tesla: Gauß-Sigma oder Grenze der Gleichverteilung.
% Die Gleichverteilung hat Standardabweichung level_t/sqrt(3).
function noise = sensor_noise(array_size, level_t, distribution)
  if strcmp(distribution, 'gaussian')
    noise = level_t * randn(array_size); % Mittelwert 0, Standardabweichung level_t.
  else
    noise = level_t * (2*rand(array_size)-1); % Grenzen +/- level_t.
  end
end


% Biot-Savart in SI-Einheiten: dB = (mu0/4pi) I (dl x r)/|r|^3.
% Die Stromlinie verlaeuft durch cfg.model_origin_m entlang der Stromrichtung.
% Jedes Segment wird durch seinen Mittelpunkt angenähert; 1e-7 ist mu0/4pi,
% 1e-6 wandelt µA in A um. Ausgabe Bx/By/Bz in Tesla.
% Berechnet das Magnetfeld eines diskretisierten geraden Stromfadens.
% ix_ua/iy_ua: Stromkomponenten / µA; x_obs/y_obs: Sensorkoordinaten / m.
% cfg liefert Geometrie und Segmentzahl; Ausgaben sind Sensormatrizen / T.
% Freiraummodell ohne Rückstrompfad; Abstand zu Segmentmitten muss > 0 sein.
function [bx, by, bz] = current_line_field(ix_ua, iy_ua, x_obs, y_obs, cfg)
  origin = [0, 0, 0];
  if isfield(cfg, 'model_origin_m'), origin = cfg.model_origin_m; end
  bx = zeros(size(x_obs)); by = bx; bz = bx;
  amplitude_ua = hypot(ix_ua, iy_ua);
  if amplitude_ua == 0, return; end
  % Richtung n=I/|I|: normierter Stromvektor. Nullstrom wurde oben abgefangen.
  nx = ix_ua / amplitude_ua; ny = iy_ua / amplitude_ua;
  % Startpunkt -L*n/2; Se5gmentvektor dl=L*n/N. Ende liegt bei +L*n/2.
  x1 = -0.5*cfg.heart_diam_m*nx; y1 = -0.5*cfg.heart_diam_m*ny;
  dlx = cfg.heart_diam_m*nx/cfg.n_segments;
  dly = cfg.heart_diam_m*ny/cfg.n_segments;
  for segment = 1:cfg.n_segments
    % Mittelpunkt des Segments; r zeigt von dieser Quelle zum Sensor.
    xs = x1+(segment-0.5)*dlx; ys = y1+(segment-0.5)*dly;
    rx = x_obs-xs-origin(1); ry = y_obs-ys-origin(2); rz = cfg.sensor_z_m-origin(3);
    r_norm = sqrt(rx.^2+ry.^2+rz.^2);
    % Kreuzprodukt für dl=(dlx,dly,0), r=(rx,ry,rz):
    % dl x r = (dly*rz, -dlx*rz, dlx*ry-dly*rx).
    % Die Beiträge aller Segmente werden nach dem Superpositionsprinzip addiert.
    factor = 1e-7 * amplitude_ua*1e-6 ./ r_norm.^3;
    bx = bx + factor .* (dly*rz);
    by = by + factor .* (-dlx*rz);
    bz = bz + factor .* (dlx.*ry-dly.*rx);
  end
end

% Richtung und Staerke ausschliesslich aus gemessenen Bx/By/Bz-Mustern.
% Kein ECG, keine wahren Stroeme und kein entfernter Offset als Fit-Eingabe.
% Inverse Musteranpassung für jede Phase: minimiert ||b-a*t(theta)||².
% sensor_values_t: Sensor-y x Sensor-x x Phase x Komponente / T.
% x_grid/y_grid: Sensorkoordinaten / m; cfg: identische Modellgeometrie.
% fit liefert Stromkomponenten / µA, Modellfeldsummen / T und RMS-Rest / T.
% Grobsuche plus lokale Winkelfeinsuche; keine Garantie eines globalen
% Minimums bei beliebigen Geometrien. Alle Messwerte haben gleiches Gewicht.
function fit = reconstruct_spatial_field(sensor_values_t, x_grid, y_grid, cfg)
  [x_obs, y_obs] = meshgrid(x_grid, y_grid);
  n_frames = size(sensor_values_t, 3);
  fit.current_ua = zeros(n_frames, 2);
  fit.magnetic_equivalent_sum_t = zeros(n_frames, 3);
  fit.residual_rms_t = zeros(n_frames, 1);
  fit.aborted = false;
  % Grobsuche: normierte Vorlagen vergleichen nur räumliche Muster.
  % 180 Grad reichen, weil die Amplitude positiv oder negativ sein darf.
  angle_step = pi/cfg.fit_angle_count;
  angles = (0:cfg.fit_angle_count-1)*angle_step;
  % Cache über mehrere Läufe: Vorlagen nur bei geänderter Geometrie,
  % Zentrierung oder Winkelraster neu berechnen. Erster Lauf kostet mehr.
  persistent cached_key cached_templates;
  origin = [0, 0, 0];
  if isfield(cfg, 'model_origin_m'), origin = cfg.model_origin_m; end
  key = [origin, x_grid(:).', y_grid(:).', cfg.heart_diam_m, cfg.sensor_z_m, ...
         cfg.n_segments, double(cfg.remove_sensor_dc), cfg.fit_angle_count];
  if isempty(cached_key) || ~isequal(cached_key, key)
    cached_templates = zeros(numel(x_obs)*3, numel(angles));
    cached_key = [];
    for k = 1:numel(angles)
      if isfield(cfg, 'should_stop')
        drawnow;
        if cfg.should_stop(), fit.aborted = true; return; end
      end
      template = spatial_model_template(angles(k), x_obs, y_obs, cfg);
      template_norm = norm(template);
      if template_norm > 0, cached_templates(:, k) = template/template_norm; end
    end
    cached_key = key;
  end
  % Jede Phase wird zu b=[Bx(:);By(:);Bz(:)] gestapelt. Die Reihenfolge
  % stimmt mit der Vorlage überein; jede Spalte ist eine Zyklusphase.
  targets = zeros(numel(x_obs)*3, n_frames);
  for k = 1:n_frames
    frame = reshape(sensor_values_t(:, :, k, :), [], 3);
    targets(:, k) = frame(:);
  end
  % Matrixprodukt berechnet alle Skalarprodukte zugleich. Bei normierten
  % Vorlagen entspricht ihr Betrag der Stärke der passenden Projektion.
  % abs berücksichtigt positive und negative Stromrichtung.
  scores = cached_templates.' * targets;
  [unused, seeds] = max(abs(scores), [], 1);
  options = optimset('TolX', cfg.fit_angle_tolerance_rad, 'Display', 'off');
  for k = 1:n_frames
    if isfield(cfg, 'should_stop')
      drawnow;
      if cfg.should_stop(), fit.aborted = true; return; end
    end
    target = targets(:, k);
    target_norm = norm(target);
    if target_norm == 0, continue; end
    % Normierung entfernt die Messamplitude aus der Winkelsuche.
    % Der Originalvektor target bleibt für die Amplitudenschätzung erhalten.
    target_direction = target/target_norm;
    seed = angles(seeds(k));
    % Feinsuche im Umfeld der besten Grobrichtung. Minimiert den nach
    % optimaler Amplitudenwahl verbleibenden quadratischen Musterfehler.
    theta = fminbnd(@(angle) spatial_fit_cost(angle, target_direction, x_obs, y_obs, cfg), ...
                    seed-angle_step, seed+angle_step, options);
    template = spatial_model_template(theta, x_obs, y_obs, cfg);
    % E=t^T*t ist die quadrierte Vorlagennorm. Nach optimaler Amplitude
  % gilt min_a ||b-a*t||² = ||b||²-(t^T*b)²/E. Der konstante erste
  % Term ist für die Winkelwahl irrelevant und wird nicht berechnet.
  energy = dot(template, template);
    if energy == 0, continue; end
    % Analytische Kleinste-Quadrate-Amplitude für festgehaltene Richtung:
    % a = (t^T b)/(t^T t), mit Vorlage t für 1 µA und Messvektor b.
    amplitude_ua = dot(template, target)/energy;
    fit.current_ua(k, :) = amplitude_ua * [cos(theta), sin(theta)];
    [bx, by, bz] = current_line_field(fit.current_ua(k, 1), fit.current_ua(k, 2), x_obs, y_obs, cfg);
    % Modellbasierte aequivalente Summe VOR raeumlicher Zentrierung.
    % Die tatsaechliche gemessene Sensorsumme bleibt bei DC-Korrektur Null.
    fit.magnetic_equivalent_sum_t(k, :) = [sum(bx(:)), sum(by(:)), sum(bz(:))];
    % RMS des unerklärten Feldes über alle Sensoren und Komponenten / T.
    % Misst Rauschen und Modellabweichung, keinen Fehler gegenüber Sollströmen.
    fit.residual_rms_t(k) = sqrt(mean((target-amplitude_ua*template).^2));
  end
  fit.method = 'least_squares_spatial_Bx_By_Bz_line_model';
  fit.equivalent_sum_definition = 'uncentered modeled sensor sum of fitted current';
end

% Dieselbe räumliche Mittelwertkorrektur wie bei den Messdaten anwenden.
% Alle Sensoren und Komponenten werden zu einem gemeinsamen Vektor gestapelt.
% Erzeugt den gestapelten Sensorfeldvektor für einen Strom von 1 µA.
% theta: Richtung / rad; Vorlage enthält Bx, dann By, dann Bz / T.
% Optional wird jede Komponente räumlich zentriert wie die Messdaten.
function template = spatial_model_template(theta, x_obs, y_obs, cfg)
  [bx, by, bz] = current_line_field(cos(theta), sin(theta), x_obs, y_obs, cfg);
  if cfg.remove_sensor_dc
    bx = bx-mean(bx(:)); by = by-mean(by(:)); bz = bz-mean(bz(:));
  end
  template = [bx(:); by(:); bz(:)];
end

% Negatives Quadrat der Projektion: Minimierung maximiert erklärtes Muster.
% Quadrieren erlaubt beide Stromvorzeichen; gleiche Gewichtung aller Sensoren.
% Zielfunktion der Winkelsuche für den normierten Messvektor.
% Eliminiert die Amplitude analytisch; Rückgabe ist die negative erklärte
% Projektionsenergie. Ein kleinerer Wert bedeutet eine bessere Musterpassung.
function cost = spatial_fit_cost(theta, target_direction, x_obs, y_obs, cfg)
  template = spatial_model_template(theta, x_obs, y_obs, cfg);
  energy = dot(template, template);
  if energy == 0, cost = 0;
  else, cost = -(dot(template, target_direction)^2)/energy; end
end

% Zeitlicher gleitender Mittelwert je Sensor und Magnetfeldkomponente.
% values: Sensor-y x Sensor-x x Phase x Komponente; Ausgabe gleicher Größe.
% window_samples: 1, 3 oder 9. Bei 1 wird die Eingabe unverändert zurückgegeben.
% Für N=2*h+1: Ausgabe(k) = sum(Eingabe(k+j), j=-h...h)/N.
% Zyklische Randbehandlung: vor Phase 1 liegt 360, nach 360 folgt 1.
% Jedes Zyklusarray wird separat gefiltert, ohne Mischung verschiedener Zyklen.
% Zentriertes Fenster vermeidet eine einseitige zeitliche Verschiebung,
% glättet aber auch schnelle Nutzsignaländerungen. Zyklische Randbehandlung
% setzt voraus, dass der Übergang zwischen Zyklusende und -anfang passend ist.
% Bei unabhängigem Rauschen sinkt dessen Standardabweichung um sqrt(N);
% benachbarte Ausgabewerte sind wegen überlappender Fenster korreliert.
function filtered = temporal_sensor_average(values, window_samples)
  if ~isscalar(window_samples) || ~ismember(window_samples, [1, 3, 9])
    error('Zeitmittel: Fensterlänge muss 1, 3 oder 9 sein.');
  end
  if window_samples == 1
    filtered = values;
    return;
  end
  half_width = (window_samples-1)/2;
  filtered = zeros(size(values));
  for shift = -half_width:half_width
    % Nur die dritte Dimension verschieben: getrennte Sensorzeitreihen.
    filtered = filtered + circshift(values, shift, 3);
  end
  filtered = filtered/window_samples;
end

% Räumlicher Modellfilter für einen Stromfaden mit genau einem Segment.
% values: Sensor-y x Sensor-x x Phase x B-Komponente / T.
% x_obs/y_obs: Sensororte / m; Quellort und Sensorhöhe stammen aus dem Modell.
% Für n_segments=1 liegt die Segmentmitte unabhängig von der Richtung am
% Modellursprung. Daher ist das Feld exakt linear in den Komponenten Ix und Iy.
% H enthält die Feldmuster für 1 µA entlang x bzw. y, jeweils Bx/By/Bz
% gestapelt. Bei aktiver DC-Korrektur sind auch die Vorlagen zentriert.
% I_hat=H\b minimiert ||b-H*I||². b_hat=H*I_hat ist die orthogonale
% Projektion auf die zwei Modellmuster: passende Signalanteile bleiben,
% orthogonales Rauschen wird entfernt. Rauschen im Modellraum bleibt bestehen.
% Keine Sollströme oder Sollspannungen als Filtereingabe; nur feste Geometrie.
% Eine anschließende Anpassung desselben Modells gewinnt dadurch keine neue
% Information. Ein kleiner Fitrest nach Projektion ist kein Genauigkeitsnachweis.
% Bei mehreren Segmenten ist dieses lineare Zwei-Muster-Modell nicht gültig.
% info enthält Stromkoeffizienten / µA und RMS des entfernten Anteils / T;
% dieser Anteil kann neben Rauschen auch reale Modellabweichungen enthalten.
function [filtered, info] = spatial_sensor_projection(values, x_obs, y_obs, cfg)
  if cfg.n_segments ~= 1
    error('Räumliche lineare Projektion verlangt n_segments=1.');
  end
  h_x = spatial_model_template(0, x_obs, y_obs, cfg);
  h_y = spatial_model_template(pi/2, x_obs, y_obs, cfg);
  model = [h_x, h_y];
  if rank(model) < 2
    error('Sensoranordnung erlaubt keine unabhängige x/y-Strombestimmung.');
  end
  dimensions = size(values);
  n_phases = size(values, 3);
  % Gleiche Stapelreihenfolge wie in der inversen Rechnung, alle Phasen
  % als Spalten. Komponenten dürfen nicht mit der Phasendimension vertauscht sein.
  frames = reshape(values, numel(x_obs), n_phases, 3);
  targets = reshape(permute(frames, [1, 3, 2]), numel(x_obs)*3, n_phases);
  coefficients = model\targets; % Stabiler Least-Squares-Fit ohne explizite Inverse.
  projected = model*coefficients;
  filtered = reshape(permute(reshape(projected, numel(x_obs), 3, n_phases), ...
      [1, 3, 2]), dimensions);
  info.enabled = true;
  info.current_ua = coefficients.';
  info.removed_rms_t = sqrt(mean((targets-projected).^2, 1)).';
  info.input_rms_t = sqrt(mean(targets.^2, 1)).';
  info.output_rms_t = sqrt(mean(projected.^2, 1)).';
  info.method = 'orthogonal_projection_on_xy_current_patterns';
end

% Leeres Array fuer Grafikhandles. MATLAB braucht Grafikobjekte (gobjects):
% numerische Handles fuehren dort z.B. bei imagesc(ax,...), title(ax,...)
% oder linkaxes zu Fehlern. Aeltere Octave-Versionen kennen gobjects nicht.
function h = empty_graphics(varargin)
  if exist('gobjects') ~= 0 %#ok<EXIST>
    h = gobjects(varargin{:});
  else
    h = zeros(varargin{:});
  end
end