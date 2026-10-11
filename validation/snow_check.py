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
WET_MM_PER_H = 0.3
SNOW_MAX_C = 1.0
WINDOW_H = 72


def _hourly(station: str, year: int) -> dict[datetime, dict]:
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


def check(acquisition: datetime) -> dict:
    """Snow/rain verdict for one acquisition time (UTC, naive or aware)."""
    t0 = acquisition.replace(tzinfo=None, minute=0, second=0, microsecond=0)
    precip = {s: _hourly(s, t0.year) for s in PRECIP_STATIONS}
    temps = {s: _hourly(s, t0.year) for s in HIGH_STATIONS}
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
                       + ", ".join(f"{HIGH_STATIONS[s]} {v:+.1f} C" for s, v in cold.items()) + ")")
    if at_acquisition:
        reasons.append("precipitating at the acquisition hour")
    totals = {s: round(sum(_value(precip[s], t, "rre150h0") or 0.0 for t in hours), 1)
              for s in precip}
    return {"passed": not reasons, "reasons": reasons, "wet_hours": len(wet_hours),
            "snow_hours": len(snow_hours), "precip_72h_mm": totals}
