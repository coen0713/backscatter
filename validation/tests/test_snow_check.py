"""Snow check: the rule on synthetic station data, and the GeoSphere reader (no network)."""

from datetime import datetime, timedelta

import pytest

from validation import pipeline, regions, snow_check


def write_geosphere(cache, station, rows):
    path = cache / "geosphere" / f"klima-v2-1h_{station}_2024.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["time,station,rr,tl,substation"]
    lines += [f"{t:%Y-%m-%dT%H:%M}+00:00,{station},{rr},{tl}," for t, rr, tl in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def tyrol_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(snow_check, "CACHE", tmp_path)
    return tmp_path


def fill(cache, acquisition, wet_hour=None, cold_hour=None):
    """Dry and warm everywhere, except an optional wet hour and cold hour."""
    top = acquisition.replace(minute=0, second=0, microsecond=0)
    hours = [top - timedelta(hours=h) for h in range(-2, 80)]
    precip, high = snow_check.STATIONS["tyrol"]
    for s in precip:
        write_geosphere(cache, s, [(t, 1.0 if t == wet_hour else 0.0, 12.0) for t in hours])
    for s in high:
        write_geosphere(cache, s, [(t, "", -2.0 if t == cold_hour else 6.0) for t in hours])


def test_geosphere_rows_map_to_the_meteoswiss_columns(tyrol_cache):
    t = datetime(2024, 8, 10, 5)
    write_geosphere(tyrol_cache, "66", [(t, 0.4, 7.5), (t + timedelta(hours=1), "", "")])
    rows = snow_check._hourly("66", 2024, "geosphere")
    assert snow_check._value(rows, t, "rre150h0") == 0.4
    assert snow_check._value(rows, t, "tre200h0") == 7.5
    assert snow_check._value(rows, t + timedelta(hours=1), "rre150h0") is None


def test_dry_scene_passes(tyrol_cache):
    acq = datetime(2024, 8, 10, 5, 30)
    fill(tyrol_cache, acq)
    r = snow_check.check(acq, "tyrol")
    assert r["passed"] and r["wet_hours"] == 0


def test_snowfall_in_the_window_fails(tyrol_cache):
    acq = datetime(2024, 8, 10, 5, 30)
    t = datetime(2024, 8, 8, 12)
    fill(tyrol_cache, acq, wet_hour=t, cold_hour=t)
    r = snow_check.check(acq, "tyrol")
    assert not r["passed"] and r["snow_hours"] == 1
    assert "Sonnblick" in r["reasons"][0]


def test_rain_at_acquisition_fails_but_warm_rain_earlier_passes(tyrol_cache):
    acq = datetime(2024, 8, 10, 5, 30)
    fill(tyrol_cache, acq, wet_hour=datetime(2024, 8, 10, 5))
    assert not snow_check.check(acq, "tyrol")["passed"]
    fill(tyrol_cache, acq, wet_hour=datetime(2024, 8, 9, 5))
    r = snow_check.check(acq, "tyrol")
    assert r["passed"] and r["wet_hours"] == 1


def test_tyrol_tiles_are_guarded():
    tiles = regions.region_tiles("tyrol")
    assert len(tiles) == 200 and tiles[0].bounds == (11.0, 46.5, 11.1, 46.6)
    if not pipeline.prereg_registered(pipeline.PREREGS["tyrol"]):
        with pytest.raises(PermissionError):
            pipeline.check_allowed(tiles[0])
