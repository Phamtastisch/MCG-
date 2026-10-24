% Hauptfunktion: steuert Datenimport, Modellbildung und Darstellung.
% Hauptfunktion ohne Argumente: liest das konfigurierte WFDB-EKG ein.
% Korrigiert Goldberger-Ableitungen, erkennt Peaks und öffnet Zeit-/Polarplots.
% Lokale Unterfunktionen übernehmen Import, Filterung, Projektion und Export.
function MCG_polar_v02()
%MCG_POLAR_V02 Einmalige Umrechnung und zwei statische Fenster.
% Lauffaehig unter MATLAB (ab R2016b, empfohlen R2020a+ wegen UTF-8-Dateien,
% Signal Processing Toolbox fuer den Hochpass) und weiterhin unter GNU Octave.
% Zeitstreifen sowie Polarplots mit Raumprojektion und CSV/MAT-Export.
% Summen je Sample aus den nullsummenkorrigierten Ableitungen, ohne Wurzel.
  % -----------------------------------------------------------------------
  % Praeambel: zentrale Parameter
  % -----------------------------------------------------------------------

  % EKG-Aufzeichnung aus der PTB-XL-Datenbank (WFDB-Format).
  % https://www.nature.com/articles/s41597-020-0495-6.
  % Nummer des EKG ..
   % cfg.rec_name          = '00003_hr'; % _hr für Samplerate 500Hz
    cfg.rec_name          = '00004_hr'; % usw.
  % heraus suchen in  https://physionet.org/content/ptb-xl/1.0.3/
  % Die Daten liegen auf der Webseite im Verzeichnis records500
  % und dann unterverzeichnis 00000 usw.
  % Weitere Angaben zur Auswahl eine EKG in ptbxl_database.csv.
  % für andere Unterverzeichnis diese Url hinten ändern, statt 00000';
   cfg.base_url          = 'https://physionet.org/files/ptb-xl/1.0.3/records500/00000';
  cfg.download_missing  = true;  % automatisches Download wenn nicht lokal schon verfuegbar
  cfg.csv_file          = 'ecg_goldberger_avr_avl_avf.csv';
  cfg.write_csv         = true;

  % Optional Gleichanteile im EKG begradigen.
  cfg.use_highpass      =  true;
  cfg.highpass_fc_hz    = 5;    %0.5;  % Cut off
  cfg.highpass_order    = 4;    % Filterordnung Butter
  cfg.peak_half_width   = 21; % x-21 ... x+21
  cfg.peak_sigma        = 7;  % Mexican-Hat-Breite in Samples.
  cfg.peak_threshold_sigma = 4; % Schwelle: Median + Faktor * robuste Streuung.
  cfg.peak_min_height_fraction = 0.5; % Mindestens halbes Maximum der Summenkurve.
  cfg.peak_pair_tolerance_fraction = 0.25; % Anteil des medianen Peakabstands.
  cfg.cycle_phase_count = 360; % Dateiformat für MCG_Kreis_v04 verlangt genau 360.
  cfg.cycle_rotation_deg = [0, 0, 0]; % Aktive Drehung: zuerst x, dann y, dann z.
  cfg.cycle_figure_pos = [100, 50, 1200, 1000];
  cfg.cycle_polar_limit_mv = 0.5; % Sichtbarer Spannungsbereich +/- Grenzwert.
  cfg.cycle_angle_grid_deg = 30; % Winkelabstand der Polar-Gitterlinien.

  % Darstellung im ersten Fenster.
  cfg.figure_id = 1;
  cfg.figure_pos = [50, 50, 1920, 1080];
  cfg.graphics_toolkit = 'qt';

  % Ab hier werden Parameter nur noch ausgewertet.
  % Mexican Hat: negative Seitenkeulen unterdrücken breite Signalanteile;
  % die positive Mitte betont schmale Peaks. Mittelwertabzug macht den
  % Filter gleichanteilsfrei; die Normierung setzt seine Energie auf eins.
  peak_x = -cfg.peak_half_width:cfg.peak_half_width;
  peak_u = peak_x / cfg.peak_sigma;
  cfg.peak_filter = (1-peak_u.^2) .* exp(-peak_u.^2/2);
  cfg.peak_filter = cfg.peak_filter - mean(cfg.peak_filter);
  cfg.peak_filter = cfg.peak_filter / norm(cfg.peak_filter);

  close all;
  clc;

  % Packages und Grafik
  setup_runtime(cfg);

  % Auslesen der Ableitungen und der parameter
  [sig_mv, fs_hz, lead_names] = read_ptb_record(cfg);
  [goldberger_mv, lead_labels] = select_leads(sig_mv, lead_names, {'AVR', 'AVL', 'AVF'});

  % Hochpass gegen Basisliniendrift; verändert auch langsame Signalanteile.
  if cfg.use_highpass
    goldberger_mv = highpass_and_zero_sum(goldberger_mv, fs_hz, cfg);
  else
    goldberger_mv = zero_sum_goldberger(goldberger_mv);
  end

  % Vektor fuer die passende Zeitachse
  n_samples = size(goldberger_mv, 1);
  % Sample 1 liegt bei t=0; jeder weitere Sampleabstand ist 1/fs_hz Sekunden.
  time_s = (0:n_samples-1).' / fs_hz;

  % CSV Datei erzeugen ...
  if cfg.write_csv
    write_goldberger_csv(cfg.csv_file, time_s, goldberger_mv);
  end

  % Quadratische Summe: stets nichtnegative Signalenergie ohne Wurzel.
  % Kubische Summe: vorzeichenabhängige Beiträge, anschließend Betrag.
  % Beide Kurven liefern getrennte Peakfolgen zur gegenseitigen Bestätigung.
  square_sum_mv2 = sum(goldberger_mv.^2, 2);
  cube_sum_mv3 = abs(sum(goldberger_mv.^3, 2)); % Betrag der kubischen Summe.
  % Vollstaendige Fenster; separate adaptive Schwelle fuer mV^2 und mV^3.
  peak_mask = false(n_samples, 2);
  power_sums = [square_sum_mv2, cube_sum_mv3];
  peak_thresholds = NaN(1, 2);
  if n_samples >= numel(cfg.peak_filter)+2
    for k = 1:2
      % 'valid' verwendet nur vollständige Filterfenster. Die Indizes
      % werden später um die halbe Filterbreite auf Samples zurückgeführt.
      response = conv(power_sums(:, k), cfg.peak_filter(:), 'valid');
      center = median(response);
      % MAD als robuste Streuung; Faktor 1.4826 entspricht bei normal-
      % verteiltem Rauschen näherungsweise der Standardabweichung.
      robust_sigma = 1.4826 * median(abs(response-center));
      % Bei konstanter Basislinie kann die MAD Null sein.
      if robust_sigma == 0
        robust_sigma = eps(max(1, max(abs(response))));
      end
      peak_thresholds(k) = center + cfg.peak_threshold_sigma*robust_sigma;
      % Lokales Maximum: links strikt größer, rechts größer oder gleich.
      % Bei einem Plateau wird so dessen linker Rand bevorzugt. +1 gleicht
      % die verkürzte Betrachtung der inneren Filterantwort aus.
      candidates = find(response(2:end-1) > peak_thresholds(k) & ...
          response(2:end-1) > response(1:end-2) & ...
          response(2:end-1) >= response(3:end)) + 1;
      % Hoehenbedingung auf der Summenkurve, nicht auf der Filterantwort.
      min_height = cfg.peak_min_height_fraction * max(power_sums(:, k));
      candidates = candidates(power_sums(candidates+cfg.peak_half_width, k) >= min_height);
      % Innerhalb einer Filter-Halbbreite nur die staerkste Antwort markieren.
      [unused, order] = sort(response(candidates), 'descend');
      selected = zeros(0, 1);
      for j = order(:).'
        candidate = candidates(j);
        if isempty(selected) || all(abs(selected-candidate) > cfg.peak_half_width)
          selected(end+1, 1) = candidate;
        end
      end
      peak_mask(selected+cfg.peak_half_width, k) = true;
      fprintf('Potenz %d: Mexican-Hat-Schwelle %.6g, %d Peaks\n', ...
          k+1, peak_thresholds(k), numel(selected));
    end
  end
  a = find(peak_mask(:, 1));
  b = find(peak_mask(:, 2));
  pairs = zeros(0, 2);
  if numel(a) >= 2 && numel(b) >= 2
    % Nur gegenseitig nächste Peaks werden gepaart; die Toleranz
    % verhindert Zuordnungen über einen zu großen Teil des Herzschlags.
    tolerance = cfg.peak_pair_tolerance_fraction*min(median(diff(a)), median(diff(b)));
    for j = 1:numel(a)
      [distance, other] = min(abs(b-a(j)));
      [unused, back] = min(abs(a-b(other)));
      if back == j && distance <= tolerance
        pairs(end+1, :) = [a(j), b(other)];
      end
    end
  end
  % Grenzen liegen zwischen bestätigten Peaks, nicht auf dem QRS-Peak.
  % Mittelung beider Peakfolgen erlaubt auch gebrochene Samplepositionen.
  mids = (pairs(1:end-1, :)+pairs(2:end, :))/2;
  boundaries = mean(mids, 2);
  % Perioden sind Abstände aufeinanderfolgender gemittelter Grenzen
  % in Samples; wegen der Mittelung nicht zwingend ganze Zahlen.
  periods = mean(diff(mids, 1, 1), 2);
  polar_fig = [];
  % Neuere MATLAB-Versionen erlauben bei figure(n) keine Parameterpaare:
  % erst Fenster mit Nummer erzeugen, dann Eigenschaften per set setzen.
  fig = figure(cfg.figure_id);
  set(fig, 'units', 'pixels', 'position', cfg.figure_pos, ...
      'name', 'MCG polar v02 - Zeitstreifen', 'numbertitle', 'off', ...
      'closerequestfcn', @quit_application);
  set(fig, 'menubar', 'none', 'toolbar', 'none');
  button_bar = uipanel('parent', fig, 'units', 'normalized', 'position', [0, 0.94, 1, 0.06]);
  uicontrol(button_bar, 'style', 'pushbutton', 'string', 'Polarplots', ...
      'units', 'normalized', 'position', [0.75, 0.15, 0.12, 0.7], ...
      'callback', @show_polar_plots);
  uicontrol(button_bar, 'style', 'pushbutton', 'string', 'Quit', ...
      'units', 'normalized', 'position', [0.89, 0.15, 0.09, 0.7], ...
      'callback', @quit_application);
  render_traces();
  show_polar_plots();

  % GUI-Callback: setzt das Abbruchflag und schließt die zugehörigen Fenster.
  % src und evt sind die vom Grafiksystem übergebenen Callback-Argumente.
  function quit_application(src, evt)
    if ~isempty(polar_fig) && ishandle(polar_fig), delete(polar_fig); end
    if ishandle(fig), delete(fig); end
  end

  % GUI-Callback: erzeugt bei Bedarf das Fenster der normierten Herzzyklen.
  % Ein vorhandenes Fenster wird aktiviert; die Zeitdaten werden nicht neu gefiltert.
  function show_polar_plots(src, evt)
    if isempty(polar_fig) || ~ishandle(polar_fig)
      polar_fig = MCG_polar_cycles(goldberger_mv, peak_mask, lead_labels, cfg);
    end
    figure(polar_fig);
    drawnow;
  end

  % Zeichnet drei Goldberger-Zeitreihen und zwei Potenzsummen.
  % Peakmarkierungen und Zwischenlinien veranschaulichen die Zyklusgrenzen;
  % die fünf Zeitachsen werden für gemeinsames Zoomen miteinander gekoppelt.
  function render_traces()
    set(0, 'currentfigure', fig);
    delete(findall(fig, 'type', 'axes'));
    traces = [goldberger_mv, square_sum_mv2, cube_sum_mv3];
    labels = {'aVR / mV', 'aVL / mV', 'aVF / mV', ...
              'aVR^2+aVL^2+aVF^2 / mV^2', '|aVR^3+aVL^3+aVF^3| / mV^3'};
    ax = empty_graphics(5, 1); % Grafikobjekte; double-Handles stoeren linkaxes in MATLAB.
    for k = 1:5
      ax(k) = axes('parent', fig, 'position', [0.12, 0.77-(k-1)*0.17, 0.84, 0.13]);
      plot(time_s, traces(:, k), 'b-');
      hold on;
      if k >= 4
        peaks = find(peak_mask(:, k-3));
        plot(time_s(peaks), traces(peaks, k), 'm*', 'markersize', 8, ...
            'linewidth', 1.2);
      end
      scale = max(abs(traces(:, k)));
      if scale == 0, scale = 1; end
      limits = 1.05*scale*[-1, 1];
      if k >= 4, limits(1) = 0; end
      hold off;
      ylim(limits);
      x_limits = [time_s(1), time_s(end)];
      if x_limits(1) == x_limits(2), x_limits = x_limits+[-0.001, 0.001]; end
      xlim(x_limits);
      grid on;
      ylabel(labels{k});
      if k == 1, title('PTB-XL: Goldberger-Ableitungen und Potenzsummen'); end
      if k == 5, xlabel('t / s'); end
    end
    % Mittelpunkte getrennt fuer beide Sternfolgen, auch halbe Samples.
    mid_samples = cell(1, 2);
    for series = 1:2
      stars = find(peak_mask(:, series));
      mid_samples{series} = (stars(1:end-1)+stars(2:end))/2;
    end
    all_mid_samples = unique([mid_samples{1}; mid_samples{2}]);
    % Gebrochene Samplepositionen mit derselben Zeitkonvention umrechnen.
    mid_times = (all_mid_samples-1)/fs_hz;
    for k = 1:5
      limits = get(ax(k), 'ylim');
      for j = 1:numel(mid_times)
        line('parent', ax(k), 'xdata', [mid_times(j), mid_times(j)], ...
            'ydata', limits, 'color', [0.65, 0.65, 0.65], ...
            'linewidth', 1, 'handlevisibility', 'off');
      end
      if k >= 4
        mids = mid_samples{k-3};
        for j = 1:numel(mids)-1
          label_time = ((mids(j)+mids(j+1))/2-1)/fs_hz;
          text(label_time, limits(2)-0.04*diff(limits), ...
              sprintf('%g Samples', mids(j+1)-mids(j)), ...
              'parent', ax(k), 'horizontalalignment', 'center', ...
              'verticalalignment', 'top', 'color', [0.35, 0.35, 0.35], ...
              'fontsize', 9, 'backgroundcolor', [1, 1, 1], 'clipping', 'on');
        end
      end
    end
    limits = get(ax(5), 'ylim');
    for cycle = 1:numel(periods)
      label_time = ((boundaries(cycle)+boundaries(cycle+1))/2-1)/fs_hz;
      text(label_time, limits(1)+0.04*diff(limits), ...
          sprintf('Mittel: %d Samples', round(periods(cycle))), ...
          'parent', ax(5), 'horizontalalignment', 'center', ...
          'verticalalignment', 'bottom', 'fontsize', 9, ...
          'backgroundcolor', [1, 1, 1], 'clipping', 'on');
    end
    linkaxes(ax, 'x');
  end
end
% Bereitet die Laufzeitumgebung vor; lädt unter Octave das Signalpaket
% für den Hochpass und wählt das konfigurierte Grafik-Toolkit.
function setup_runtime(cfg)
  if ~is_octave()
    % MATLAB: butter und filtfilt stammen aus der Signal Processing Toolbox.
    if cfg.use_highpass && (isempty(which('butter')) || isempty(which('filtfilt')))
      error(['Signal Processing Toolbox (butter/filtfilt) nicht gefunden. ', ...
          'Toolbox installieren oder cfg.use_highpass = false setzen.']);
    end
    return;
  end
  if is_octave()
    if cfg.use_highpass
    try
      eval('pkg load signal');
    catch
      error('Octave-Paket "signal" wird benoetigt: pkg install -forge signal');
    end

    end

    if ~isempty(cfg.graphics_toolkit)
      try
        graphics_toolkit(cfg.graphics_toolkit);
      catch
        warning('Grafik-Toolkit "%s" konnte nicht gesetzt werden.', cfg.graphics_toolkit);
      end
    end
  end
end

% Prueft, ob der Code unter GNU Octave laeuft.
% Liefert true, wenn die eingebaute OCTAVE_VERSION-Konstante existiert.
% Erlaubt Octave-spezifische Vorbereitung bei gemeinsamem MATLAB-Code.
function tf = is_octave()
  tf = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
end

% Liest eine PTB-XL/WFDB-Aufzeichnung und skaliert die Rohdaten nach mV.
% Liest Header und Binärdaten der konfigurierten WFDB-Aufzeichnung.
% Ausgaben: Samples x Kanäle / mV, Abtastrate / Hz und Kanalnamen.
function [sig_mv, fs_hz, lead_names] = read_ptb_record(cfg)
  ensure_record_files(cfg);

  hea_file = [cfg.rec_name '.hea'];
  dat_file = [cfg.rec_name '.dat'];

  [n_chan, fs_hz, n_samples, gain, baseline, lead_names] = parse_wfdb_header(hea_file);
  raw = read_wfdb_dat(dat_file, n_chan, n_samples);

  % WFDB-Format: Messwert = (Rohwert - Baseline) / Gain, hier in mV.
  sig_mv = bsxfun(@rdivide, bsxfun(@minus, raw, baseline), gain);
end

% Stellt sicher, dass Header- und Datendatei lokal vorhanden sind.
% Prüft, ob .hea und .dat lokal vorhanden sind.
% Lädt fehlende Dateien bei erlaubtem Download aus cfg.base_url nach.
function ensure_record_files(cfg)
  for ext_cell = {'hea', 'dat'}
    ext = ext_cell{1};
    file_name = [cfg.rec_name '.' ext];

    if exist(file_name, 'file')
      continue;
    end

    if ~cfg.download_missing
      error('Datei fehlt: %s', file_name);
    end

    url = [cfg.base_url '/' file_name];
    % Plattformunabhaengiger Download ohne externes wget:
    % MATLAB nutzt websave, Octave urlwrite.
    try
      if is_octave()
        urlwrite(url, file_name);
      else
        websave(file_name, url, weboptions('Timeout', 60));
      end
    catch err
      if exist(file_name, 'file'), delete(file_name); end
      error('Download fehlgeschlagen: %s (%s)', url, err.message);
    end
  end
end

% Extrahiert Kanalzahl, Abtastrate, Kalibrierung und Kanalnamen aus dem WFDB-Header.
% Liest Kanalzahl, Abtastrate / Hz, Samplezahl und Kanalnamen.
% Gain ist die Skalierung in ADC-Einheiten/mV; baseline der ADC-Nullwert.
% Erwartet das hier verwendete WFDB-Kalibrierformat gain(baseline)/mV.
function [n_chan, fs_hz, n_samples, gain, baseline, lead_names] = parse_wfdb_header(hea_file)
  lines = regexp(fileread(hea_file), '\r\n|\n|\r', 'split');
  first = strsplit(strtrim(lines{1}));

  n_chan   = str2double(first{2});
  fs_hz    = str2double(first{3});
  n_samples = str2double(first{4});

  gain = zeros(1, n_chan);
  baseline = zeros(1, n_chan);
  lead_names = cell(1, n_chan);

  for k = 1:n_chan
    parts = strsplit(strtrim(lines{k+1}));
    cal_text = parts{3};
    tok = regexp(cal_text, '([+-]?[0-9.]+)\((-?[0-9]+)\)/mV', 'tokens');

    if isempty(tok)
      error('Kalibrierung in Headerzeile %d nicht erkannt: %s', k+1, cal_text);
    end

    gain(k) = str2double(tok{1}{1});
    baseline(k) = str2double(tok{1}{2});
    lead_names{k} = parts{end};
  end
end

% Liest die gespeicherten int16-Samples aus der EKG-Datendatei (_wfdb-Datei).
% Liest verschachtelte 16-Bit-Kanäle mit Little-Endian-Bytefolge.
% Rückgabe: Samples x Kanäle als double, zunächst noch ADC-Rohwerte.
% Gilt für das hier verwendete Format 16, nicht für sämtliche WFDB-Formate.
function raw = read_wfdb_dat(dat_file, n_chan, n_samples)
  fid = fopen(dat_file, 'rb');
  if fid < 0
    error('Kann Datei nicht oeffnen: %s', dat_file);
  end

  cleaner = onCleanup(@() fclose(fid));
  % fread füllt spaltenweise: eine Spalte enthält alle Kanäle eines Samples.
  % Transponieren liefert die übliche Matrix Samples x Kanäle.
  raw = fread(fid, [n_chan, n_samples], 'int16=>double', 0, 'ieee-le');
  raw = raw.';

  if size(raw, 1) ~= n_samples || size(raw, 2) ~= n_chan
    error('Unerwartete Datengroesse in %s.', dat_file);
  end
end

% Waehlt die benoetigten Ableitungen in der vorgegebenen Reihenfolge aus.
% Wählt Kanäle nach Namen in der gewünschten Reihenfolge aus.
% sig_mv: Samples x Kanäle; sel_mv enthält nur die gesuchten Ableitungen.
% Ein fehlender Kanal führt zu einem Fehler statt einer stillen Ersatzwahl.
function [sel_mv, labels] = select_leads(sig_mv, lead_names, labels)
  idx = zeros(1, numel(labels));

  for k = 1:numel(labels)
    idx_k = find(strcmp(lead_names, labels{k}), 1);
    if isempty(idx_k)
      error('Ableitung fehlt im Datensatz: %s', labels{k});
    end
    idx(k) = idx_k;
  end

  sel_mv = sig_mv(:, idx);
end

% Entfernt langsame Basisliniendrift und erzwingt anschliessend die Nullsumme der Goldberger-Ableitungen
% Entfernt langsame Drift mittels Butterworth-Hochpass und filtfilt.
% Eingabe/Ausgabe: Samples x Ableitungen / mV; fs_hz: Abtastrate / Hz.
% Anschließend wird die Goldberger-Nullsummenbedingung wiederhergestellt.
function sig_filt = highpass_and_zero_sum(sig_mv, fs_hz, cfg)
  % butter erwartet die auf die Nyquistfrequenz fs/2 normierte Grenze.
  % Nur 0 < wn < 1 ist als digitale Hochpassgrenze zulässig.
  wn = cfg.highpass_fc_hz / (fs_hz / 2);
  if wn <= 0 || wn >= 1
    error('Ungueltige Hochpass-Grenzfrequenz: %.3g Hz', cfg.highpass_fc_hz);
  end

  % Butterworth: glatter Durchlassbereich. Vorwärts-/Rückwärtsfilterung
  % mit filtfilt vermeidet Phasenverschiebung; die effektive Ordnung
  % verdoppelt sich und an den Signalrändern bleiben Filtereffekte möglich.
  [b, a] = butter(cfg.highpass_order, wn, 'high');
  sig_filt = filtfilt(b, a, sig_mv);
  sig_filt = zero_sum_goldberger(sig_filt);
end

% Korrigiert die Goldberger-Ableitungen auf eine zeilenweise Summe von Null
% gleichzeitige EKG-Sample sollen zusammen Null werden.
% Zieht je Sample den Mittelwert der drei Ableitungen von allen drei ab.
% Dies ist die kleinste quadratische Korrektur zur Bedingung Summe=0.
% Eingabe/Ausgabe: Samples x drei Ableitungen in derselben Einheit.
function sig_corr = zero_sum_goldberger(sig_mv)
  % Fuer ideale Goldberger-Ableitungen gilt aVR + aVL + aVF = 0.
  % Orthogonale Projektion auf die Ebene aVR+aVL+aVF=0:
  % gleicher Korrekturwert für alle Kanäle, Unterschiede bleiben erhalten.
  sig_corr = sig_mv - mean(sig_mv, 2);
end

% Schreibt Zeitachse und Goldberger-Ableitungen als CSV-Datei.
% Schreibt je Sample Zeit / s und aVR/aVL/aVF / mV in eine CSV-Datei.
% onCleanup schließt die Datei auch bei einem Schreibfehler.
function write_goldberger_csv(csv_file, time_s, goldberger_mv)
  fid = fopen(csv_file, 'w');
  if fid < 0
    error('Kann CSV-Datei nicht schreiben: %s', csv_file);
  end

  cleaner = onCleanup(@() fclose(fid));
  fprintf(fid, 't_s,aVR_mV,aVL_mV,aVF_mV\n');
  % fprintf verarbeitet spaltenweise. out transponieren, damit jeweils
  % Zeit und drei Kanäle desselben Samples als eine CSV-Zeile erscheinen.
  out = [time_s, goldberger_mv];
  fprintf(fid, '%.9g,%.9g,%.9g,%.9g\n', out.');
end



% Lokale Unterfunktion: verbindet die beiden Peakfolgen und interpoliert
% vollständige Herzzyklen auf ein gemeinsames Phasenraster.
% goldberger_mv: Samples x drei Ableitungen / mV; peak_mask: zwei logische
% Peakspalten; lead_labels: Kanalnamen; cfg: Vorgaben der Hauptfunktion.
% Rückgabe fig ist das Fensterhandle für Projektion, Winkelwahl und Export.
function fig = MCG_polar_cycles(goldberger_mv, peak_mask, lead_labels, cfg)
% Lokale Unterfunktion: Zyklusnormierung, Raumprojektion und Export.
% Verwendet ausschließlich die am Anfang der Hauptfunktion gesetzten Parameter.
  rotation_deg = cfg.cycle_rotation_deg;
  a = find(peak_mask(:, 1));
  b = find(peak_mask(:, 2));
  pairs = zeros(0, 2);
  if numel(a) >= 2 && numel(b) >= 2
    tolerance = cfg.peak_pair_tolerance_fraction*min(median(diff(a)), median(diff(b)));
    for j = 1:numel(a)
      [distance, other] = min(abs(b-a(j)));
      [unused, back] = min(abs(a-b(other)));
      if back == j && distance <= tolerance
        pairs(end+1, :) = [a(j), b(other)];
      end
    end
  end
  mids = (pairs(1:end-1, :)+pairs(2:end, :))/2;
  boundaries = mean(mids, 2);
  periods = mean(diff(mids, 1, 1), 2);
  fig = figure('name', 'Goldberger: Polarplots der Herzzyklen', ...
      'numbertitle', 'off', 'position', cfg.cycle_figure_pos, 'visible', 'on');
  if isempty(periods)
    axis off;
    text(0.5, 0.5, 'Zu wenige passende Peaks fuer einen vollstaendigen Herzzyklus.', ...
        'horizontalalignment', 'center');
    drawnow;
    return;
  end
  % Linear auf gleiche Phasen interpolieren: Start und Ende sind enthalten.
  % Unterschiedliche Zyklusdauern werden so vergleichbar; die Winkelachse
  % beschreibt die Zyklusphase, nicht die Richtung des Herzvektors.
  phase_deg = linspace(1, 360, cfg.cycle_phase_count);
  values = zeros(numel(periods), cfg.cycle_phase_count, 3);
  for cycle = 1:numel(periods)
    % phase=1 entspricht der linken Grenze, phase=360 der rechten.
    % (phase-1)/359 läuft von 0 bis 1; Multiplikation mit der Periodenlänge
    % liefert die Sampleposition. interp1 gewichtet die zwei Nachbarsamples
    % linear; die Originaldaten werden dabei zeitlich gestaucht oder gedehnt.
    query = boundaries(cycle)+(phase_deg-1)/359*periods(cycle);
    for lead = 1:3
      values(cycle, :, lead) = interp1((1:size(goldberger_mv, 1)).', ...
          goldberger_mv(:, lead), query, 'linear');
    end
  end
  if numel(rotation_deg) ~= 3 || any(~isfinite(rotation_deg))
    delete(fig);
    error('Drei endliche Drehwinkel in Grad erforderlich.');
  end
  projected_xyz = zeros(size(values));
  angle_bar = uipanel('parent', fig, 'units', 'normalized', ...
      'position', [0, 0.92, 1, 0.08]);
  angle_inputs = empty_graphics(1, 3);
  axis_names = {'x', 'y', 'z'};
  for dimension = 1:3
    left = 0.01+(dimension-1)*0.19;
    uicontrol('parent', angle_bar, 'style', 'text', ...
        'string', sprintf('Winkel %s (Grad)', axis_names{dimension}), ...
        'units', 'normalized', 'position', [left, 0.24, 0.11, 0.45], ...
        'horizontalalignment', 'left');
    angle_inputs(dimension) = uicontrol('parent', angle_bar, 'style', 'edit', ...
        'string', num2str(rotation_deg(dimension)), 'tag', ['rotation_' axis_names{dimension}], ...
        'units', 'normalized', 'position', [left+0.11, 0.20, 0.065, 0.55], ...
        'callback', @change_rotation);
  end
  uicontrol('parent', angle_bar, 'style', 'pushbutton', 'string', 'Anwenden', ...
      'units', 'normalized', 'position', [0.59, 0.20, 0.11, 0.55], ...
      'callback', @change_rotation);
  uicontrol('parent', angle_bar, 'style', 'pushbutton', 'string', 'CSV speichern', ...
      'tag', 'save_csv', 'units', 'normalized', 'position', [0.72, 0.20, 0.13, 0.55], ...
      'callback', @(src, evt) save_xyz('csv'));
  uicontrol('parent', angle_bar, 'style', 'pushbutton', 'string', 'MAT speichern', ...
      'tag', 'save_mat', 'units', 'normalized', 'position', [0.86, 0.20, 0.13, 0.55], ...
      'callback', @(src, evt) save_xyz('mat'));
  render_projection();

  % GUI-Callback: exportiert die zuletzt gedrehten, unbeschnittenen x/y/z-Daten.
  % format bestimmt CSV oder MAT; MAT speichert zusätzlich Zyklusmetadaten.
  % CSV enthält eine Zeile pro Zyklusphase, MAT eine Zeile pro Zyklus.
  function save_xyz(format)
    [filename, directory] = uiputfile(['*.' format], 'x/y/z-Herzzyklen speichern', ...
        ['MCG_xyz_' datestr(now, 'yyyymmdd_HHMMSS') '.' format]);
    if isequal(filename, 0), return; end
    [unused, unused_name, extension] = fileparts(filename);
    if isempty(extension), filename = [filename '.' format]; end
    target = fullfile(directory, filename);
    % Zeilen: Zyklen; Spalten: 360 Phasenwerte. Unbeschnittene Spannungen.
    x_mv = projected_xyz(:, :, 1);
    y_mv = projected_xyz(:, :, 2);
    z_mv = projected_xyz(:, :, 3);
    mean_xyz_mv = reshape(mean(projected_xyz, 1), cfg.cycle_phase_count, 3);
    period_samples = periods;
    cycle_boundary_samples = boundaries;
    try
      if strcmp(format, 'mat')
        save(target, 'x_mv', 'y_mv', 'z_mv', 'mean_xyz_mv', 'phase_deg', ...
            'period_samples', 'cycle_boundary_samples', 'rotation_deg', '-v7');
      else
        fid = fopen(target, 'w');
        if fid < 0, error('Die Datei konnte nicht geoeffnet werden.'); end
        try
          fprintf(fid, 'cycle,phase_deg,x_mV,y_mV,z_mV,rotation_x_deg,rotation_y_deg,rotation_z_deg\n');
          for cycle = 1:numel(periods)
            data = [repmat(cycle, cfg.cycle_phase_count, 1), phase_deg(:), ...
                reshape(projected_xyz(cycle, :, :), cfg.cycle_phase_count, 3), ...
                repmat(rotation_deg(:).', cfg.cycle_phase_count, 1)];
            fprintf(fid, '%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n', data.');
          end
        catch write_err
          fclose(fid);
          rethrow(write_err);
        end
        fclose(fid); % Datei vor Rueckkehr aus dem GUI-Callback schliessen.
      end
      fprintf('x/y/z gespeichert: %s\n', target);
    catch err
      errordlg(err.message, 'Speichern fehlgeschlagen');
    end
  end

  % GUI-Callback: liest drei endliche Winkel / Grad aus den Eingabefeldern.
  % Aktualisiert die Raumprojektion und damit die anschließend exportierten Daten.
  function change_rotation(src, evt)
    new_angles = zeros(1, 3);
    for dimension = 1:3
      new_angles(dimension) = str2double(get(angle_inputs(dimension), 'string'));
    end
    if any(~isfinite(new_angles))
      errordlg('Bitte drei endliche Zahlen eingeben.', 'Ungueltige Winkel');
      return;
    end
    rotation_deg = new_angles;
    render_projection();
  end

  % Rechnet Goldberger-Ableitungen in einen ebenen x/y-Vektor um.
  % Dreht diesen im Raum und zeichnet Ableitungen sowie x/y/z über Zyklusphase.
  % Die Mittelkurve entsteht durch Mittelung an gleicher normierter Phase.
  function render_projection()
  set(0, 'currentfigure', fig);
  delete(findall(fig, 'type', 'axes'));
  % Goldberger-Achsen: aVR=(-sqrt(3)/2,+1/2),
  % aVL=(+sqrt(3)/2,+1/2), aVF=(0,-1). Ausgangs-z ist Null.
  xyz = zeros(size(values));
  % Inversion der drei Projektionen bei Nullsumme:
  % aVL-aVR=sqrt(3)*x, aVR+aVL=y und aVF=-y. Die symmetrische
  % y-Formel nutzt alle drei Ableitungen. z=0, da nur die Frontalebene vorliegt.
  xyz(:, :, 1) = (values(:, :, 2)-values(:, :, 1))/sqrt(3);
  xyz(:, :, 2) = (values(:, :, 1)+values(:, :, 2)-2*values(:, :, 3))/3;
  % Trigonometrische Funktionen erwarten Radiant: Grad * pi/180.
  angles = rotation_deg*pi/180;
  cx = cos(angles(1)); sx = sin(angles(1));
  cy = cos(angles(2)); sy = sin(angles(2));
  cz = cos(angles(3)); sz = sin(angles(3));
  rx = [1,0,0; 0,cx,-sx; 0,sx,cx];
  ry = [cy,0,sy; 0,1,0; -sy,0,cy];
  rz = [cz,-sz,0; sz,cz,0; 0,0,1];
  % Rechtshändige aktive Drehung: x, dann y, dann z. Das Transponieren
  % ist nötig, weil die Vektoren hier als Zeilen gespeichert sind.
  % z entsteht durch Drehung der Frontalebene, nicht durch neue Messdaten.
  rotation = rz*ry*rx;
  % Alle Zyklus-/Phasenvektoren zu Zeilen zusammenfassen, gemeinsam drehen
  % und wieder in die Form Zyklus x Phase x Raumkomponente zurückformen.
  xyz = reshape(reshape(xyz, [], 3)*rotation.', size(xyz));
  projected_xyz = xyz;
  all_values = cat(3, values(:, :, [2, 1, 3]), xyz);
  plot_labels = [reshape(lead_labels([2, 1, 3]), 1, []), {'x', 'y', 'z'}];
  % Feste Spannungsgrenzen: Zentrum -0.5 mV, Nullkreis, Aussenring +0.5 mV.
  % Radius = Spannung + Grenzwert: negative Spannungen liegen innerhalb
  % des Nullkreises. NaN blendet Werte außerhalb des Sichtbereichs aus;
  % die exportierten Daten bleiben vollständig und unbeschnitten.
  offset = cfg.cycle_polar_limit_mv;
  theta = phase_deg*pi/180;
  for lead = 1:6
    col = mod(lead-1, 3);
    row = floor((lead-1)/3);
    axes('parent', fig, 'position', [0.045+0.325*col, 0.51-0.43*row, 0.27, 0.32]);
    v = all_values(:, :, lead);
    % Eigenes Polargitter: Beschriftungen zeigen Spannung statt Offsetradius.
    axis equal;
    axis off;
    xlim(2.45*offset*[-1, 1]);
    ylim(2.45*offset*[-1, 1]);
    hold on;
    circle_theta = linspace(0, 2*pi, 361);
    patch(offset*cos(circle_theta(:)), offset*sin(circle_theta(:)), [0.85, 0.85, 0.85], ...
        'edgecolor', 'none', 'handlevisibility', 'off');
    voltage_ticks = offset*[-0.5, 0, 0.5, 1];
    for voltage = voltage_ticks
      radius = voltage+offset;
      if voltage == 0
        plot(radius*cos(circle_theta), radius*sin(circle_theta), 'k--', 'linewidth', 1.2);
      else
        plot(radius*cos(circle_theta), radius*sin(circle_theta), '-', 'color', [0.65, 0.65, 0.65]);
      end
    end
    for grid_angle = 0:cfg.cycle_angle_grid_deg:(360-cfg.cycle_angle_grid_deg)
      phi = grid_angle*pi/180;
      plot([0, 2*offset*cos(phi)], [0, 2*offset*sin(phi)], '-', 'color', [0.7, 0.7, 0.7]);
      degree_label = grid_angle;
      if grid_angle == 0, degree_label = 360; end
      text(2.2*offset*cos(phi), 2.2*offset*sin(phi), sprintf('%d°', degree_label), ...
          'horizontalalignment', 'center', 'fontsize', 9);
    end
    for cycle = 1:size(v, 1)
      radius = v(cycle, :)+offset;
      radius(radius < 0 | radius > 2*offset) = NaN;
      plot(radius.*cos(theta), radius.*sin(theta), 'b-');
    end
    % Phasenmittel statt zeitlichem Mittel: gleiche Phase verschiedener
    % Zyklen mitteln, anschließend denselben Spannungs-Offset addieren.
    radius = mean(v, 1)+offset;
    radius(radius < 0 | radius > 2*offset) = NaN;
    h_mean = plot(radius.*cos(theta), radius.*sin(theta), 'r-');
    set(h_mean, 'linewidth', 2);
    % Spannungslabels zuletzt zeichnen, damit Kurven sie nicht verdecken.
    for voltage = voltage_ticks
      radius = voltage+offset;
      text(radius*cos(pi/4), radius*sin(pi/4), sprintf('%+.3g mV', voltage), ...
          'horizontalalignment', 'center', 'fontsize', 9, 'backgroundcolor', [1, 1, 1]);
    end
    text(0, -2.65*offset, 'Hellgrau: negativ | ausserhalb 0 mV: positiv', ...
        'horizontalalignment', 'center', 'fontsize', 9);
    hold off;
    title(sprintf('%s: %d Herzzyklen', plot_labels{lead}, numel(periods)));
  end
  set(fig, 'name', sprintf('Goldberger und xyz | Drehung x=%g, y=%g, z=%g Grad', rotation_deg));
  drawnow;
  end
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