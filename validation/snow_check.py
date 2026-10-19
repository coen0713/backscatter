"""Snow (and rain) check for candidate scenes, from MeteoSwiss open data.

docs/prereg.md, section 2: a scene is rejected if the MeteoSwiss record shows
fresh snow above 2,500 m within the 3 days before it. Operational definition
(fixed before any imagery was obtained), on hourly station data:

* Precipitation evidence: >= 0.3 mm in an hour at Guetsch (GUE, 2286 m) or
  Pilatus (PIL, 2105 m). No MeteoSwiss gauge above 2,500 m inside the region
  reports in 2024 (Titlis has no 2024 precipitation values).
* Snow evidence: during such an hour, any station at or above 2,500 m in or
  next to the region is at or below +1.0 C: Jungfraujoch (JUN, 3571 m),
  Gornergrat (GOR, 3129 m), Diablerets (DIA, 2964 m), Eggishorn (EGH, 2892 m),
  Les Attelas (ATT, 2734 m), Titlis (TIT, 3045 m). Snow can fall at up to
  about +1 C, so this flags snowfall at that station's altitude.
* Reject if snow evidence occurs in the 72 hours before acquisition, or if it
  is precipitating at the acquisition hour (wet surfaces change C-band
  backscatter by several dB, independent of snow).

Hours are checked from acquisition - 72 h to acquisition + 1 h, so either
convention for MeteoSwiss hourly time stamps (interval start or end) is
covered. Limits: two gauges at ~2,100-2,300 m stand in for the whole region,
and stations below the highest summits (4,000+ m) cannot see snowfall there.
The SLF IMIS network would measure snow height directly above 2,500 m but its
open API only reaches 7 days back.

Data: https://data.geo.admin.ch, collection ch.meteoschweiz.ogd-smn (open
government data, no account), cached in data/cache/meteoswiss.
"""

from __future__ import annotations

import csv
from datetime import datetime, timedelta
from pathlib import Path

from validation.fetch import CACHE, download

BASE = "https://data.geo.admin.ch/ch.meteoschweiz.ogd-smn"
PRECIP_STATIONS = {"gue": "Guetsch 2286 m", "pil": "Pilatus 2105 m"}
HIGH_STATIONS = {"jun": "Jungfraujoch 3571 m", "gor": "Gornergrat 3129 m",
                 "dia": "Diablerets 2964 m", "egh": "Eggishorn 2892 m",
                 "att": "Les Attelas 2734 m", "tit": "Titlis 3045 m"}
# Station sets per region: (precipitation gauges, temperature stations >= 2,500 m).
# "bernese" is round 1 (docs/prereg.md); "graubuenden" is round 2
# (docs/prereg_angular.md, section 8), fixed before its scene selection.
STATIONS = {
    "bernese": (PRECIP_STATIONS, HIGH_STATIONS),
    "graubuenden": (
        {"sbe": "San Bernardino 1639 m", "pma": "Piz Martegnas 2668 m", "vab": "Valbella 1568 m",
         "biv": "Bivio 1856 m", "aro": "Arosa 1878 m", "sia": "Segl-Maria 1804 m",
         "wfj": "Weissfluhjoch 2691 m", "cov": "Piz Corvatsch 3294 m", "dav": "Davos 1594 m",
         "sam": "Samedan 1709 m", "beh": "Passo del Bernina 2260 m", "buf": "Buffalora 1971 m"},
        {"pma": "Piz Martegnas 2668 m", "wfj": "Weissfluhjoch 2691 m",
         "cov": "Piz Corvatsch 3294 m"}),
    # Round 3 (docs/prereg_partb.md, section 8): GeoSphere Austria station ids,
    # chosen before its scene selection. Gauges: every station at 1,500 m or
    # higher in the region with complete 2024 precipitation; temperature:
    # stations at 2,500 m or higher in or next to it.
    "tyrol": (
        {"82": "Rudolfshuette 2317 m", "196": "Patscherkofel 2251 m",
         "85": "Schmittenhoehe 1956 m", "66": "Obergurgl 1941 m",
         "12216": "Hahnenkamm-Sonnenrast 1667 m", "15430": "Kolm-Saigurn 1626 m",
         "14912": "Hintertux 1505 m"},
        {"213": "Sonnblick 3109 m", "17320": "Brunnenkogel 3437 m",
         "184": "Pitztaler Gletscher 2864 m"}),
}
# Where each region's stations come from (MeteoSwiss SwissMetNet or GeoSphere Austria).
SOURCES = {"bernese": "meteoswiss", "graubuenden": "meteoswiss", "tyrol": "geosphere"}
GEOSPHERE = "https://dataset.api.hub.geosphere.at/v1/station/historical/klima-v2-1h"
WET_MM_PER_H = 0.3
SNOW_MAX_C = 1.0
WINDOW_H = 72


def _hourly_geosphere(station: str, year: int, retries: int = 8) -> dict[datetime, dict]:
    """GeoSphere Austria hourly data (open, no account) as rows keyed like the
    MeteoSwiss ones: rr (mm) -> rre150h0, tl (deg C) -> tre200h0. Times are UTC."""
    path = CACHE / "geosphere" / f"klima-v2-1h_{station}_{year}.csv"
    if not path.exists():
        # GeoSphere's certificate chain is not in Python's default store on
        # every system, so use certifi's bundle (a rasterio dependency).
        import ssl  # noqa: PLC0415
        import time  # noqa: PLC0415
        import urllib.error  # noqa: PLC0415
        import urllib.request  # noqa: PLC0415

        import certifi  # noqa: PLC0415
        url = (f"{GEOSPHERE}?parameters=rr,tl&station_ids={station}&start={year}-01-01T00:00"
               f"&end={year + 1}-01-01T00:00&output_format=csv")
        ctx = ssl.create_default_context(cafile=certifi.where())
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(url, timeout=300, context=ctx) as resp:
                    text = resp.read().decode("utf-8")
                break
            except (urllib.error.URLError, ConnectionError, TimeoutError):
                if attempt == retries - 1:
                    raise
                time.sleep(min(2.0 ** attempt, 30.0))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    out = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            t = datetime.fromisoformat(row["time"]).replace(tzinfo=None)
            if t.year == year:
                out[t] = {"rre150h0": row["rr"], "tre200h0": row["tl"]}
    return out


def _hourly(station: str, year: int, source: str = "meteoswiss") -> dict[datetime, dict]:
    if source == "geosphere":
        return _hourly_geosphere(station, year)
    decade = year // 10 * 10
    name = f"ogd-smn_{station}_h_historical_{decade}-{decade + 9}.csv"
    path = download(f"{BASE}/{station}/{name}", CACHE / "meteoswiss" / name)
    out = {}
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f, delimiter=";"):
            t = datetime.strptime(row["reference_timestamp"], "%d.%m.%Y %H:%M")
            if t.year == year:
                out[t] = row
    return out


def _value(table: dict[datetime, dict], t: datetime, column: str) -> float | None:
    v = table.get(t, {}).get(column)
    return float(v) if v not in (None, "") else None


def check(acquisition: datetime, region: str = "bernese") -> dict:
    """Snow/rain verdict for one acquisition time (UTC, naive or aware)."""
    precip_names, high_names = STATIONS[region]
    t0 = acquisition.replace(tzinfo=None, minute=0, second=0, microsecond=0)
    source = SOURCES[region]
    precip = {s: _hourly(s, t0.year, source) for s in precip_names}
    temps = {s: _hourly(s, t0.year, source) for s in high_names}
    hours = [t0 + timedelta(hours=h) for h in range(-WINDOW_H, 2)]
    wet_hours, snow_hours = [], []
    for t in hours:
        p = max((_value(precip[s], t, "rre150h0") or 0.0) for s in precip)
        if p < WET_MM_PER_H:
            continue
        wet_hours.append(t)
        cold = {s: v for s in temps if (v := _value(temps[s], t, "tre200h0")) is not None
                and v <= SNOW_MAX_C}
        if cold:
            snow_hours.append((t, cold))
    at_acquisition = [t for t in wet_hours if t in (t0, t0 + timedelta(hours=1))]
    reasons = []
    if snow_hours:
        t, cold = snow_hours[0]
        reasons.append(f"snowfall evidence in {len(snow_hours)} h (first {t:%Y-%m-%d %H:%M}: "
                       + ", ".join(f"{high_names[s]} {v:+.1f} C" for s, v in cold.items()) + ")")
    if at_acquisition:
        reasons.append("precipitating at the acquisition hour")
    totals = {s: round(sum(_value(precip[s], t, "rre150h0") or 0.0 for t in hours), 1)
              for s in precip}
    return {"passed": not reasons, "reasons": reasons, "wet_hours": len(wet_hours),
            "snow_hours": len(snow_hours), "precip_72h_mm": totals}
