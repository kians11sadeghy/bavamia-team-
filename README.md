# Lunar South Pole Mission Navigator

A mission-planning web dashboard and MATLAB analysis tool for checking **sunlight hours**, **solar power potential**, and **Direct-to-Earth (DTE) communication windows** at lunar south pole landing sites.

Created for the **NASA International Space Apps Challenge**.

---

## Why We Built This
At the Moon's south pole, the Sun stays very close to the horizon. Even small crater rims can cast long shadows or block line-of-sight communication with Earth. When launch dates slip or mission plans change, teams have to re-check lighting and communication conditions across different sites—which usually takes a lot of manual work in tools like JPL Horizons and Quickmap.

We built this tool to make comparing landing sites and dates fast, visual, and easy to understand.

---

## What It Does
- **Preset Sites & Custom Coordinates:** Comes with a built-in list of candidate south pole landing sites (*Shackleton Rim*, *Mons Mouton*, *Malapert Massif*, *Nobile Rim*, *Haworth Rim*), and also lets you enter any custom latitude and longitude.
- **Side-by-Side Comparison:** Compare two landing sites or two different launch dates on the same screen.
- **Interactive Timeline Slider:** Scrub hour-by-hour through the mission window to watch the Sun and Earth move across the local horizon, with live countdowns to the next sunrise/sunset and Earthrise/Earthset.
- **Energy & Comms Metrics:** Calculates total sunlit time, longest continuous darkness period, a solar power generation proxy (`Power Proxy`), and Earth visibility percentage.

---

## Files in This Repo
- `index.html` — Standalone web app (runs directly in any browser with no installation) showing the sky/horizon view, visibility timeline, and comparison table.
- `lunar_site_eval.m` — MATLAB script that pulls ephemeris data from the **NASA JPL Horizons API** (`coord@301`), caches results locally, accounts for Sun (0.265 deg) and Earth (0.95 deg) angular radii, and calculates visibility and power metrics.
- `lunar_make_mask.m` — MATLAB script that reads **LRO LOLA South Polar GeoTIFF DEM** elevation files (`LDEM_80S_..._ADJ.TIF`) and accounts for lunar curvature (R = 1737.4 km) to build a 360-degree horizon profile.
- `lunar_navigator_app.m` — MATLAB GUI wrapper for running the site evaluation and terrain mask tools interactively.

---

## Data & Tools Used
- **Ephemeris Data:** [NASA JPL Horizons System](https://ssd.jpl.nasa.gov/horizons/)
- **Lunar Topography:** [NASA LRO LOLA Elevation Data (DEM)](https://pds-geosciences.wustl.edu/missions/lro/lola.htm) via Lunar Quickmap
- **Development Note:** We used **Claude AI** to help write parts of the code and build the interactive web interface alongside our MATLAB calculations.

---

## How to Run
- **Web App:** Just open `index.html` in your browser, or use the live GitHub Pages link.
- **MATLAB:** Place the `.m` files in your MATLAB path and run:
  `lunar_site_eval(-89.54, 0.0, 'Name', 'Shackleton Rim', 'Start', '2027-03-10', 'Days', 30)`
