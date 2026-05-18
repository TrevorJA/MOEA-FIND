"""Tests for the unified slug builder (src.io_paths.slugs).

Added with the Phase-7 stage-04 slug migration: stage 04 now builds
output-path slugs via ``moea_slug`` instead of the retired
``make_variant_slug``. Asserts the slug is deterministic, stable, and
round-trips through ``parse_slug``.
"""

from src.io_paths.slugs import moea_slug, parse_slug, build_slug


class TestMoeaSlug:
    def test_basic_shape_and_roundtrip(self):
        s = moea_slug(mode="residual", n_years=20, nfe=200000, seed=42,
                       ssi=3, metrics="primary", cons="dv-l2")
        assert s.startswith("moea__")
        p = parse_slug(s)
        assert p["_stage"] == "moea"
        assert p["mode"] == "residual"
        assert int(p["T"]) == 20
        assert p["nfe"] == 200000          # nfe auto-parsed back to int
        assert int(p["s"]) == 42
        assert p["ssi"] == "3"
        assert p["metrics"] == "primary"
        assert p["cons"] == "dv-l2"

    def test_deterministic(self):
        kw = dict(mode="index", n_years=10, nfe=1000, seed=7, ssi=3)
        assert moea_slug(**kw) == moea_slug(**kw)

    def test_extra_sorted_and_present(self):
        s = moea_slug(mode="residual", n_years=10, nfe=5000, seed=1,
                       cons="dv-l2", extra={"st": "ad", "sfx": "v2"})
        # extra keys appear, seed token stays last.
        assert "st=ad" in s and "sfx=v2" in s
        assert s.rsplit("__", 1)[1] == "s=1"

    def test_distinct_configs_distinct_slugs(self):
        a = moea_slug(mode="residual", n_years=10, nfe=1000, seed=1)
        b = moea_slug(mode="residual", n_years=10, nfe=1000, seed=2)
        assert a != b

    def test_none_optionals_omitted(self):
        s = moea_slug(mode="residual", n_years=10, nfe=1000, seed=1)
        assert "ssi=" not in s and "metrics=" not in s and "cons=" not in s


class TestBuildSlugContract:
    def test_stage_prefix(self):
        assert build_slug("moea", mode="x", s=1).startswith("moea__")
