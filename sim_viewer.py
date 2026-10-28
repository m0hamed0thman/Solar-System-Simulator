#!/usr/bin/env python3
"""
sim_viewer.py -- Solar System Simulator viewer (reads traj.csv from planets_mpi.c)

Usage:   python3 sim_viewer.py [traj.csv] [--epoch YYYY-MM-DD]

--epoch   calendar date of t = 0 in the file. If the CSV has a comment line such as
          "# epoch=2026-01-01" it is detected automatically. Without it the viewer
          shows "Day N" instead of calendar dates.

Mouse:    hover a body for details | click a body to focus on it | wheel = zoom
Keyboard: Space play/pause | Left/Right step | + / - zoom | Home full system
          T trails | L labels | S snapshot | F fullscreen | Q quit
"""

import sys
import csv
import re
import time
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import to_rgba
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle

# We build our own toolbar/shortcuts, so switch off Matplotlib's defaults.
for _key in list(plt.rcParams.keys()):
    if _key.startswith("keymap."):
        try:
            plt.rcParams[_key] = []
        except Exception:
            pass
plt.rcParams["toolbar"] = "None"

# ---------------------------------------------------------------------------
# Constants and body information
# ---------------------------------------------------------------------------
AU_KM = 149597870.7
LIGHT_MIN_PER_AU = 8.3167
SECONDS_PER_DAY = 86400.0

# Colors (dark professional theme)
BG = "#0a0e17"
SPACE = "#04060c"
PANEL = "#141b2b"
BTN = "#1a2338"
BTN_HOVER = "#2a3a5c"
EDGE = "#2b3650"
TEXT = "#d5dcea"
MUTED = "#7f8aa3"
ACCENT = "#4da3ff"

RADII_KM = {
    'Sun': 696000, 'Mercury': 2439.7, 'Venus': 6051.8, 'Earth': 6371.0,
    'Moon': 1737.4, 'Mars': 3389.5, 'Phobos': 11.27, 'Deimos': 6.2,
    'Jupiter': 69911, 'Io': 1821.6, 'Europa': 1560.8,
    'Ganymede': 2634.1, 'Callisto': 2410.3,
    'Saturn': 58232, 'Mimas': 198.2, 'Enceladus': 252.1,
    'Tethys': 531.1, 'Dione': 561.4, 'Rhea': 763.8,
    'Titan': 2574.7, 'Iapetus': 734.5,
    'Uranus': 25362, 'Ariel': 578.9, 'Umbriel': 584.7,
    'Titania': 788.4, 'Oberon': 761.4, 'Miranda': 235.8,
    'Neptune': 24622, 'Triton': 1353.4,
    'Pluto': 1188.3, 'Charon': 606.0,
    'Ceres': 473.0, 'Pallas': 256.0, 'Vesta': 262.7,
    'Eros': 8.42, 'Psyche': 113.0, 'Davida': 144.7, 'Interamnia': 166.0,
}

PARENT = {
    'Moon': 'Earth', 'Phobos': 'Mars', 'Deimos': 'Mars',
    'Io': 'Jupiter', 'Europa': 'Jupiter', 'Ganymede': 'Jupiter', 'Callisto': 'Jupiter',
    'Mimas': 'Saturn', 'Enceladus': 'Saturn', 'Tethys': 'Saturn', 'Dione': 'Saturn',
    'Rhea': 'Saturn', 'Titan': 'Saturn', 'Iapetus': 'Saturn',
    'Ariel': 'Uranus', 'Umbriel': 'Uranus', 'Titania': 'Uranus',
    'Oberon': 'Uranus', 'Miranda': 'Uranus',
    'Triton': 'Neptune', 'Charon': 'Pluto',
}

PLANETS = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
DWARFS = ['Pluto', 'Ceres']

COLORS = {
    'Sun': '#FDB813', 'Mercury': '#B1ADAD', 'Venus': '#E8CDA2', 'Earth': '#4F8FDF',
    'Moon': '#C8C8C8', 'Mars': '#D1603D', 'Phobos': '#A08070', 'Deimos': '#B09080',
    'Jupiter': '#D9A066', 'Io': '#E6D36F', 'Europa': '#D8C8A8',
    'Ganymede': '#A39A90', 'Callisto': '#7A6E62',
    'Saturn': '#E3CF9A', 'Mimas': '#BBBBBB', 'Enceladus': '#F2F2F2',
    'Tethys': '#D9D9D9', 'Dione': '#CFCFCF', 'Rhea': '#C0C0C0',
    'Titan': '#E0A030', 'Iapetus': '#A09580',
    'Uranus': '#8AD4E0', 'Ariel': '#B8C8D0', 'Umbriel': '#8A98A0',
    'Titania': '#B0BCC4', 'Oberon': '#98A4AC', 'Miranda': '#C4D0D8',
    'Neptune': '#4766E8', 'Triton': '#D8BFD8',
    'Pluto': '#C8A890', 'Charon': '#9A9A9A',
    'Ceres': '#B8A89A', 'Pallas': '#A89888', 'Vesta': '#C0B098',
    'Eros': '#B09078', 'Psyche': '#A8A8B0', 'Davida': '#988878', 'Interamnia': '#A09080',
}

# Quick views: (button text, body to center on, view radius in AU)
VIEWS = [
    ("Full system",    'Sun',     None),     # radius computed from the data
    ("Inner planets",  'Sun',     1.8),
    ("Asteroid belt",  'Sun',     4.0),
    ("Gas giants",     'Sun',     11.0),
    ("Earth & Moon",   'Earth',   0.004),
    ("Mars system",    'Mars',    0.0002),
    ("Jupiter system", 'Jupiter', 0.015),
    ("Saturn system",  'Saturn',  0.028),
    ("Uranus system",  'Uranus',  0.0048),
    ("Neptune system", 'Neptune', 0.003),
    ("Pluto & Charon", 'Pluto',   0.00016),
]

TRAIL_MODES = [("off", 0), ("30 d", 30), ("1 yr", 365), ("full", None)]


# ---------------------------------------------------------------------------
# Small helper functions
# ---------------------------------------------------------------------------
def body_kind(name):
    if name == 'Sun':
        return "Star"
    if name in PLANETS:
        return "Planet"
    if name in DWARFS:
        return "Dwarf planet"
    if name in PARENT:
        return "Moon of " + PARENT[name]
    return "Asteroid"


def visual_points(radius_km):
    """Marker diameter in points: log-scaled so every body is visible
    and bigger bodies always look bigger than smaller ones."""
    x = np.log10(max(radius_km, 5.0) / 5.0) / np.log10(696000.0 / 5.0)
    return 2.0 + 12.0 * x ** 2.2


def nice_number(x):
    """Round x to 1, 2 or 5 times a power of ten."""
    exponent = np.floor(np.log10(x))
    f = x / 10 ** exponent
    if f < 1.5:
        n = 1
    elif f < 3.5:
        n = 2
    elif f < 7.5:
        n = 5
    else:
        n = 10
    return n * 10 ** exponent


def nice_length_au(span_au):
    """A round length close to span_au (round in AU for big spans, in km for small)."""
    if span_au >= 0.1:
        return nice_number(span_au)
    return nice_number(span_au * AU_KM) / AU_KM


def format_length(au):
    km = au * AU_KM
    if au >= 0.1:
        return "%g AU" % au
    if km >= 1e6:
        return "%g million km" % (km / 1e6)
    return "{:,.0f} km".format(km)


def format_speed(days_per_second):
    if days_per_second < 1.0:
        return "%.1f h/s" % (days_per_second * 24.0)
    if days_per_second < 365.0:
        return "%.1f d/s" % days_per_second
    return "%.1f yr/s" % (days_per_second / 365.25)


# ---------------------------------------------------------------------------
# Loading the data
# ---------------------------------------------------------------------------
def load_traj(filename):
    """Returns names, times (N,), positions (N, bodies, 3), epoch text or None."""
    header = None
    epoch_text = None
    body_names, body_index = [], {}
    t_list, b_list, xyz_list = [], [], []

    with open(filename) as f:
        for line in csv.reader(f):
            if not line:
                continue
            if line[0].startswith('#'):
                joined = ",".join(line)
                match = re.search(r"epoch\D*(\d{4}-\d{2}-\d{2})", joined)
                if match and epoch_text is None:
                    epoch_text = match.group(1)
                continue
            if header is None:
                header = [h.strip() for h in line]
                continue
            row = dict(zip(header, line))
            name = row['name'].strip()
            if name not in body_index:
                body_index[name] = len(body_names)
                body_names.append(name)
            t_list.append(float(row['t_days']))
            b_list.append(body_index[name])
            xyz_list.append((float(row['x_au']), float(row['y_au']), float(row['z_au'])))

    if not t_list:
        raise SystemExit(filename + ": no data rows")

    times = np.unique(t_list)
    ti = np.searchsorted(times, t_list)
    positions = np.full((len(times), len(body_names), 3), np.nan)
    positions[ti, b_list] = xyz_list
    return body_names, times, positions, epoch_text


# ---------------------------------------------------------------------------
# The viewer
# ---------------------------------------------------------------------------
class Viewer:
    def __init__(self, names, times, pos, epoch):
        self.names = names
        self.B = len(names)
        self.times = times
        self.pos = pos
        self.epoch = epoch                      # datetime or None
        self.index = {n: i for i, n in enumerate(names)}

        self.radius_km = np.array([RADII_KM.get(n, 100.0) for n in names])
        self.radius_au = self.radius_km / AU_KM
        self.visual_pts = np.array([visual_points(r) for r in self.radius_km])
        self.colors = [COLORS.get(n, '#9aa5b8') for n in names]
        self.draw_order = np.argsort(-self.radius_km)          # big first, small on top
        self.priority = list(np.argsort(-self.radius_km))      # label priority
        self.label_fs = [9.0 if (n == 'Sun' or n in PLANETS) else 7.5 for n in names]

        self.sun = self.index.get('Sun', 0)
        self.center = self.sun
        self.full_radius = self.auto_full_radius()
        self.radius = self.full_radius

        self.t = float(times[0])
        self.playing = True
        self.speed = 1.0                        # days per second
        self.labels_on = True
        self.trail_mode = 1
        self.active_view = 0
        self.guard = False                      # stops slider callbacks looping
        self.toast_until = 0.0
        self.last_clock = time.perf_counter()
        self.last_rel = np.zeros((self.B, 3))
        self.last_frame = 0
        self.widgets = []

        self.build_figure()
        self.build_controls()
        self.connect_events()
        self.draw()

    # ---------- data helpers ----------
    def auto_full_radius(self):
        finite = self.pos[np.isfinite(self.pos)]
        max_r = np.max(np.abs(finite)) if finite.size else 35.0
        return float(min(80.0, max(35.0, max_r * 1.15)))

    def positions_at(self, t):
        """Linear interpolation between the two nearest saved frames."""
        i = int(np.searchsorted(self.times, t, side='right')) - 1
        i = max(0, min(i, len(self.times) - 2))
        span = self.times[i + 1] - self.times[i]
        f = (t - self.times[i]) / span
        f = min(max(f, 0.0), 1.0)
        return self.pos[i] * (1 - f) + self.pos[i + 1] * f, i

    def date_text(self):
        if self.epoch is None:
            return "Day %.2f" % self.t
        when = self.epoch + timedelta(days=self.t)
        return when.strftime("%Y-%m-%d  %H:%M")

    # ---------- figure ----------
    def build_figure(self):
        self.fig = plt.figure(figsize=(14, 8.4), facecolor=BG)
        try:
            self.fig.canvas.manager.set_window_title("Solar System Simulator")
        except Exception:
            pass

        self.ax = self.fig.add_axes([0.14, 0.19, 0.85, 0.78], facecolor=SPACE)
        ax = self.ax
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(EDGE)

        # background stars (fixed in the window, not in space)
        rng = np.random.default_rng(7)
        ax.scatter(rng.random(350), rng.random(350), s=rng.random(350) * 1.6 + 0.2,
                   c='white', alpha=0.35, transform=ax.transAxes, zorder=0, linewidths=0)

        # distance rings (centered on the focus body)
        self.rings, self.ring_texts = [], []
        for _ in range(8):
            ring = Circle((0, 0), 1, fill=False, ec='#2a3652', lw=0.7,
                          ls=(0, (3, 4)), zorder=1)
            ax.add_patch(ring)
            self.rings.append(ring)
            self.ring_texts.append(ax.text(0, 0, "", color='#55627f', fontsize=7,
                                           zorder=2, ha='left', va='bottom'))

        # trails, Sun glow, bodies
        self.trail_lines = LineCollection([], linewidths=0.8, zorder=5)
        ax.add_collection(self.trail_lines)
        self.trail_lines.set_color([to_rgba(c, 0.55) for c in self.colors])

        self.glow = ax.scatter([0], [0], s=1, c=COLORS['Sun'], alpha=0.16,
                               linewidths=0, zorder=6)
        order = self.draw_order
        self.scat = ax.scatter(np.zeros(self.B), np.zeros(self.B),
                               s=np.ones(self.B), c=[self.colors[i] for i in order],
                               edgecolors='#000000', linewidths=0.4, zorder=10)

        # labels sit directly above each body (offset in points, so it works at any zoom)
        outline = [pe.withStroke(linewidth=2.4, foreground=SPACE)]
        self.labels = []
        for b, name in enumerate(self.names):
            ann = ax.annotate(name, xy=(0, 0), xytext=(0, 6), textcoords='offset points',
                              ha='center', va='bottom', fontsize=self.label_fs[b],
                              color=self.colors[b], zorder=30, path_effects=outline)
            self.labels.append(ann)

        # info panel, scale bar, toast, hover tooltip
        self.info = ax.text(0.012, 0.985, "", transform=ax.transAxes, va='top',
                            family='monospace', fontsize=9, color=TEXT, zorder=50,
                            bbox=dict(boxstyle='round,pad=0.6', fc=BG, ec=EDGE, alpha=0.82))
        self.bar, = ax.plot([0, 1], [0, 0], color=TEXT, lw=1.8, zorder=40,
                            solid_capstyle='butt')
        self.bar_text = ax.text(0, 0, "", color=TEXT, fontsize=8, ha='left',
                                va='bottom', zorder=40)
        self.toast = ax.text(0.5, 0.04, "", transform=ax.transAxes, ha='center',
                             color='white', fontsize=10, zorder=60,
                             bbox=dict(boxstyle='round,pad=0.5', fc=ACCENT, ec='none'))
        self.toast.set_visible(False)
        self.tip = ax.annotate("", xy=(0, 0), xytext=(14, 14), textcoords='offset points',
                               fontsize=8.5, color=TEXT, family='monospace', zorder=70,
                               bbox=dict(boxstyle='round,pad=0.5', fc=PANEL, ec=ACCENT, alpha=0.96))
        self.tip.set_visible(False)

    # ---------- controls ----------
    def add_button(self, rect, text, callback, fs=9):
        bax = self.fig.add_axes(rect)
        btn = Button(bax, text, color=BTN, hovercolor=BTN_HOVER)
        btn.label.set_color(TEXT)
        btn.label.set_fontsize(fs)
        for spine in bax.spines.values():
            spine.set_color(EDGE)
        btn.on_clicked(lambda _e: callback())
        self.widgets.append(btn)
        return btn

    def add_slider(self, rect, label, vmin, vmax, vinit, callback):
        sax = self.fig.add_axes(rect, facecolor=PANEL)
        sl = Slider(sax, label, vmin, vmax, valinit=vinit, color=ACCENT)
        sl.label.set_color(TEXT)
        sl.label.set_fontsize(9)
        sl.valtext.set_visible(False)
        for spine in sax.spines.values():
            spine.set_color(EDGE)
        sl.on_changed(callback)
        self.widgets.append(sl)
        return sl

    def build_controls(self):
        fig = self.fig

        # sidebar: quick views
        fig.text(0.07, 0.955, "VIEWS", color=ACCENT, fontsize=10, fontweight='bold',
                 ha='center', va='center')
        self.view_buttons = []
        y = 0.895
        for i, (text, _body, _r) in enumerate(VIEWS):
            btn = self.add_button([0.008, y, 0.12, 0.043], text,
                                  lambda i=i: self.set_view(i), fs=8.5)
            self.view_buttons.append(btn)
            y -= 0.052
        help_lines = ["Space  play / pause", "< >    step frame", "+ -    zoom",
                      "Wheel  zoom", "Click  focus body", "Hover  details",
                      "T trails   L labels", "S snapshot", "Home   full system"]
        fig.text(0.012, y - 0.01, "\n".join(help_lines), color=MUTED, fontsize=7.5,
                 va='top', family='monospace')

        # sliders
        self.time_slider = self.add_slider([0.20, 0.135, 0.70, 0.022], "Time",
                                           self.times[0], self.times[-1], self.t,
                                           self.on_time_slider)
        self.speed_slider = self.add_slider([0.20, 0.095, 0.30, 0.022], "Speed",
                                            -2.0, 3.0, 0.0, self.on_speed_slider)
        self.zoom_slider = self.add_slider([0.62, 0.095, 0.28, 0.022], "Zoom",
                                           np.log10(2e-5), 2.0, np.log10(self.radius),
                                           self.on_zoom_slider)
        self.time_text = fig.text(0.91, 0.146, "", color=TEXT, fontsize=8.5, va='center')
        self.speed_text = fig.text(0.505, 0.106, "", color=TEXT, fontsize=8.5, va='center')
        self.zoom_text = fig.text(0.905, 0.106, "", color=TEXT, fontsize=8.5, va='center')

        # button row
        self.play_btn = self.add_button([0.140, 0.030, 0.065, 0.040], "Pause", self.toggle_play)
        self.add_button([0.210, 0.030, 0.035, 0.040], "<", lambda: self.step(-1))
        self.add_button([0.250, 0.030, 0.035, 0.040], ">", lambda: self.step(1))
        self.trail_btn = self.add_button([0.295, 0.030, 0.095, 0.040], "", self.cycle_trails)
        self.label_btn = self.add_button([0.395, 0.030, 0.095, 0.040], "", self.toggle_labels)
        self.add_button([0.495, 0.030, 0.085, 0.040], "Snapshot", self.snapshot)

        tax = self.fig.add_axes([0.640, 0.030, 0.120, 0.040])
        self.goto_box = TextBox(tax, "Go to  ", initial="", color=PANEL, hovercolor=BTN_HOVER)
        self.goto_box.label.set_color(TEXT)
        self.goto_box.label.set_fontsize(9)
        self.goto_box.text_disp.set_color(TEXT)
        self.goto_box.on_submit(self.goto)
        self.widgets.append(self.goto_box)

        presets = [("1 h/s", 1 / 24.0), ("1 d/s", 1.0), ("10 d/s", 10.0), ("1 yr/s", 365.25)]
        x = 0.800
        for text, value in presets:
            self.add_button([x, 0.030, 0.046, 0.040], text,
                            lambda v=value: self.set_speed(v), fs=8)
            x += 0.050

        self.refresh_button_texts()
        self.highlight_view()

    def connect_events(self):
        c = self.fig.canvas
        c.mpl_connect('key_press_event', self.on_key)
        c.mpl_connect('scroll_event', self.on_scroll)
        c.mpl_connect('motion_notify_event', self.on_motion)
        c.mpl_connect('button_press_event', self.on_click)
        c.mpl_connect('resize_event', lambda _e: self.draw())
        try:
            self.anim = FuncAnimation(self.fig, self.on_tick, interval=33,
                                      cache_frame_data=False)
        except TypeError:
            self.anim = FuncAnimation(self.fig, self.on_tick, interval=33)

    # ---------- actions ----------
    def toast_message(self, text, seconds=2.0):
        self.toast.set_text(text)
        self.toast.set_visible(True)
        self.toast_until = time.perf_counter() + seconds

    def refresh_button_texts(self):
        self.play_btn.label.set_text("Pause" if self.playing else "Play")
        self.trail_btn.label.set_text("Trails: " + TRAIL_MODES[self.trail_mode][0])
        self.label_btn.label.set_text("Labels: " + ("on" if self.labels_on else "off"))

    def highlight_view(self):
        for i, btn in enumerate(self.view_buttons):
            color = ACCENT if i == self.active_view else BTN
            btn.color = color
            btn.ax.set_facecolor(color)
            btn.label.set_color('black' if i == self.active_view else TEXT)

    def toggle_play(self):
        self.playing = not self.playing
        self.last_clock = time.perf_counter()
        self.refresh_button_texts()

    def step(self, direction):
        self.playing = False
        i = int(np.searchsorted(self.times, self.t, side='right')) - 1
        if direction < 0 and self.t <= self.times[i] + 1e-9:
            i -= 1
        j = min(max(i + (1 if direction > 0 else 0), 0), len(self.times) - 1)
        self.set_time(self.times[j])
        self.refresh_button_texts()

    def set_time(self, t):
        self.t = float(min(max(t, self.times[0]), self.times[-1]))
        self.guard = True
        self.time_slider.set_val(self.t)
        self.guard = False
        self.draw()

    def set_speed(self, days_per_second):
        self.speed = days_per_second
        self.guard = True
        self.speed_slider.set_val(np.log10(days_per_second))
        self.guard = False
        self.draw()

    def set_radius(self, radius):
        self.radius = float(min(max(radius, 2e-5), 100.0))
        self.guard = True
        self.zoom_slider.set_val(np.log10(self.radius))
        self.guard = False

    def set_view(self, i):
        text, body, radius = VIEWS[i]
        if body not in self.index:
            self.toast_message(body + " is not in the data file")
            return
        self.center = self.index[body]
        self.set_radius(self.full_radius if radius is None else radius)
        self.active_view = i
        self.highlight_view()
        self.draw()

    def focus_on(self, b):
        name = self.names[b]
        for i, (_t, body, radius) in enumerate(VIEWS):
            if body == name and radius is not None and name != 'Sun':
                self.set_view(i)
                return
        if name == 'Sun':
            self.set_view(0)
            return
        self.center = b
        if name in PARENT and PARENT[name] in self.index:
            P, _ = self.positions_at(self.t)
            d = np.linalg.norm(P[b] - P[self.index[PARENT[name]]])
            self.set_radius(2.5 * d if np.isfinite(d) and d > 0 else 0.001)
        else:
            self.set_radius(0.5)
        self.active_view = None
        self.highlight_view()
        self.draw()

    def cycle_trails(self):
        self.trail_mode = (self.trail_mode + 1) % len(TRAIL_MODES)
        self.refresh_button_texts()
        self.draw()

    def toggle_labels(self):
        self.labels_on = not self.labels_on
        self.refresh_button_texts()
        self.draw()

    def snapshot(self):
        stamp = self.date_text().replace(" ", "_").replace(":", "")
        name = "snapshot_%s.png" % stamp
        extent = self.ax.get_window_extent().transformed(self.fig.dpi_scale_trans.inverted())
        self.fig.savefig(name, bbox_inches=extent, dpi=150, facecolor=BG)
        self.toast_message("Saved " + name)
        print("saved", name)

    def goto(self, text):
        text = text.strip()
        if not text:
            return
        try:
            if re.match(r"\d{4}-\d{2}-\d{2}", text) and self.epoch is not None:
                when = datetime.strptime(text[:10], "%Y-%m-%d")
                t = (when - self.epoch).total_seconds() / SECONDS_PER_DAY
            else:
                t = float(text)
        except ValueError:
            self.toast_message("Type a day number" + (" or YYYY-MM-DD" if self.epoch else ""))
            return
        if t < self.times[0] or t > self.times[-1]:
            self.toast_message("Out of range (%.1f to %.1f days)" % (self.times[0], self.times[-1]))
        self.set_time(t)

    # ---------- widget callbacks ----------
    def on_time_slider(self, value):
        if self.guard:
            return
        self.t = float(value)
        self.draw()

    def on_speed_slider(self, value):
        if self.guard:
            return
        self.speed = 10 ** value
        self.draw()

    def on_zoom_slider(self, value):
        if self.guard:
            return
        self.radius = 10 ** value
        self.active_view = None
        self.highlight_view()
        self.draw()

    # ---------- mouse and keyboard ----------
    def on_key(self, event):
        key = event.key
        if key == ' ':
            self.toggle_play()
        elif key == 'right':
            self.step(1)
        elif key == 'left':
            self.step(-1)
        elif key in ('+', '='):
            self.set_radius(self.radius / 1.25)
            self.draw()
        elif key in ('-', '_'):
            self.set_radius(self.radius * 1.25)
            self.draw()
        elif key == 'home':
            self.set_view(0)
        elif key == 't':
            self.cycle_trails()
        elif key == 'l':
            self.toggle_labels()
        elif key == 's':
            self.snapshot()
        elif key == 'f':
            try:
                self.fig.canvas.manager.full_screen_toggle()
            except Exception:
                pass
        elif key in ('q', 'escape'):
            plt.close(self.fig)

    def on_scroll(self, event):
        if event.inaxes is not self.ax:
            return
        factor = 0.8 if event.button == 'up' else 1.25
        self.set_radius(self.radius * factor)
        self.draw()

    def nearest_body(self, event):
        disp = self.ax.transData.transform(self.last_rel[:, :2])
        dist = np.hypot(disp[:, 0] - event.x, disp[:, 1] - event.y)
        dist[~np.isfinite(dist)] = np.inf
        b = int(np.argmin(dist))
        limit = max(9.0, self.last_sizes[b] / 2 * self.fig.dpi / 72 + 5)
        return b if dist[b] <= limit else None

    def on_motion(self, event):
        if event.inaxes is not self.ax:
            if self.tip.get_visible():
                self.tip.set_visible(False)
                self.fig.canvas.draw_idle()
            return
        b = self.nearest_body(event)
        if b is None:
            if self.tip.get_visible():
                self.tip.set_visible(False)
                self.fig.canvas.draw_idle()
            return
        self.tip.xy = (self.last_rel[b, 0], self.last_rel[b, 1])
        self.tip.set_text(self.tooltip_text(b))
        self.tip.set_visible(True)
        self.fig.canvas.draw_idle()

    def on_click(self, event):
        if event.inaxes is not self.ax or event.button != 1:
            return
        b = self.nearest_body(event)
        if b is not None:
            self.focus_on(b)

    def tooltip_text(self, b):
        name = self.names[b]
        i = self.last_frame
        P, _ = self.positions_at(self.t)
        lines = [name, body_kind(name)]
        if name != 'Sun':
            d = np.linalg.norm(P[b] - P[self.sun])
            lines.append("From Sun : %.4f AU (%.1f light-min)" % (d, d * LIGHT_MIN_PER_AU))
        ref = self.index.get(PARENT.get(name, 'Sun'), self.sun)
        if name in PARENT and ref != self.sun:
            d = np.linalg.norm(P[b] - P[ref]) * AU_KM
            lines.append("From %s: {:,.0f} km".format(d) % self.names[ref])
        if b != ref:
            a = (self.pos[i, b] - self.pos[i, ref])
            c = (self.pos[i + 1, b] - self.pos[i + 1, ref])
            dt = (self.times[i + 1] - self.times[i]) * SECONDS_PER_DAY
            v = np.linalg.norm(c - a) / dt * AU_KM
            if np.isfinite(v):
                lines.append("Speed vs %s: %.2f km/s" % (self.names[ref], v))
        lines.append("Radius   : {:,.1f} km".format(self.radius_km[b]))
        return "\n".join(lines)

    # ---------- animation tick ----------
    def on_tick(self, _frame):
        now = time.perf_counter()
        dt = min(now - self.last_clock, 0.1)
        self.last_clock = now

        if self.toast.get_visible() and now > self.toast_until:
            self.toast.set_visible(False)

        if self.playing:
            self.t += self.speed * dt
            if self.t > self.times[-1]:
                self.t = float(self.times[0])
            self.guard = True
            self.time_slider.set_val(self.t)
            self.guard = False
            self.draw()

    # ---------- drawing ----------
    def draw(self):
        ax = self.ax
        P, frame = self.positions_at(self.t)
        self.last_frame = frame

        center = P[self.center]
        if not np.all(np.isfinite(center)):
            center = np.zeros(3)
        rel = P - center
        self.last_rel = rel

        # view window: vertical half-size = radius, horizontal follows the window shape
        box = ax.get_window_extent()
        aspect = box.width / max(box.height, 1.0)
        R = self.radius
        ax.set_xlim(-R * aspect, R * aspect)
        ax.set_ylim(-R, R)

        # marker sizes: log-scaled look, but never smaller than the true size
        pts_per_au = box.height / (2 * R) * 72.0 / self.fig.dpi
        true_pts = 2 * self.radius_au * pts_per_au
        sizes = np.minimum(np.maximum(self.visual_pts, true_pts), 900.0)
        self.last_sizes = sizes

        order = self.draw_order
        self.scat.set_offsets(rel[order, :2])
        self.scat.set_sizes(sizes[order] ** 2)
        self.glow.set_offsets([rel[self.sun, :2]])
        self.glow.set_sizes([(sizes[self.sun] * 3.2) ** 2])

        self.draw_labels(rel, sizes, box)
        self.draw_trails(rel, frame)
        self.draw_rings(R, aspect)
        self.draw_scale_bar(R, aspect)
        self.draw_info()
        self.fig.canvas.draw_idle()

    def draw_labels(self, rel, sizes, box):
        if not self.labels_on:
            for ann in self.labels:
                ann.set_visible(False)
            return
        disp = self.ax.transData.transform(rel[:, :2])
        dpi = self.fig.dpi
        taken = []
        order = [self.center] + [b for b in self.priority if b != self.center]
        for b in order:
            ann = self.labels[b]
            x, y = disp[b]
            ok = (np.isfinite(x) and np.isfinite(y) and
                  box.x0 <= x <= box.x1 and box.y0 <= y <= box.y1)
            if ok:
                fs_px = self.label_fs[b] * dpi / 72.0
                w = 0.6 * fs_px * len(self.names[b]) + 4
                h = 1.4 * fs_px
                off_pts = sizes[b] / 2 + 3
                cx, cy = x, y + off_pts * dpi / 72.0 + h / 2
                for tx, ty, tw, th in taken:
                    if abs(cx - tx) < (w + tw) / 2 and abs(cy - ty) < (h + th) / 2:
                        ok = False
                        break
            if ok:
                taken.append((cx, cy, w, h))
                ann.xy = (rel[b, 0], rel[b, 1])
                ann.set_position((0, off_pts))
                ann.set_visible(True)
            else:
                ann.set_visible(False)

    def draw_trails(self, rel, frame):
        days = TRAIL_MODES[self.trail_mode][1]
        if days == 0:
            self.trail_lines.set_segments([])
            return
        i1 = frame + 1
        i0 = 0 if days is None else int(np.searchsorted(self.times, self.t - days))
        i0 = min(i0, frame)
        block = self.pos[i0:i1 + 1]
        seg = block - block[:, self.center:self.center + 1]
        stride = max(1, seg.shape[0] // 700)
        segments = []
        for b in range(self.B):
            line = np.vstack([seg[::stride, b, :2], rel[b:b + 1, :2]])
            segments.append(line)
        self.trail_lines.set_segments(segments)

    def draw_rings(self, R, aspect):
        step = nice_length_au(R / 2.5)
        reach = R * aspect * 1.05
        for k, (ring, text) in enumerate(zip(self.rings, self.ring_texts), start=1):
            radius = k * step
            if radius > reach:
                ring.set_visible(False)
                text.set_visible(False)
                continue
            ring.set_radius(radius)
            ring.set_visible(True)
            x, y = radius * np.cos(0.5), radius * np.sin(0.5)
            text.set_position((x, y))
            text.set_text(format_length(radius))
            text.set_visible(abs(x) < R * aspect and abs(y) < R)

    def draw_scale_bar(self, R, aspect):
        length = nice_length_au(R * 0.35)
        x0 = -R * aspect * 0.96
        y0 = -R * 0.92
        self.bar.set_data([x0, x0 + length], [y0, y0])
        self.bar_text.set_position((x0, y0 + R * 0.02))
        self.bar_text.set_text(format_length(length))

    def draw_info(self):
        lines = [
            "SOLAR SYSTEM SIMULATOR",
            self.date_text(),
            "Focus : %s" % self.names[self.center],
            "View  : +/- %s" % format_length(self.radius),
            "Speed : %s" % format_speed(self.speed),
            "Bodies: %d" % self.B,
        ]
        self.info.set_text("\n".join(lines))
        self.time_text.set_text(self.date_text())
        self.speed_text.set_text(format_speed(self.speed))
        self.zoom_text.set_text(format_length(self.radius))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = sys.argv[1:]
    epoch_arg = None
    if '--epoch' in args:
        k = args.index('--epoch')
        epoch_arg = args[k + 1]
        del args[k:k + 2]
    filename = args[0] if args else 'traj.csv'

    names, times, pos, epoch_text = load_traj(filename)
    if len(times) < 2:
        raise SystemExit("need at least 2 frames")

    epoch_text = epoch_arg or epoch_text
    epoch = datetime.strptime(epoch_text, "%Y-%m-%d") if epoch_text else None

    print("loaded %d bodies, %d frames, t = %.2f .. %.2f days"
          % (len(names), len(times), times[0], times[-1]))
    for k, name in enumerate(names):
        if not np.all(np.isfinite(pos[:, k])):
            print("WARNING: '%s' has NaN values (simulation may have exploded)" % name)

    viewer = Viewer(names, times, pos, epoch)
    plt.show()


if __name__ == '__main__':
    main()