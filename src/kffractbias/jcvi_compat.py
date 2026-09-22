"""Run JCVI catalog with quota help compatible with Python 3.14 argparse."""

from __future__ import annotations

from argparse import Action
from typing import Any


def main() -> None:
    from jcvi.apps.base import OptionParser
    from jcvi.compara import catalog, quota

    class QuotaOptionParser(OptionParser):
        def add_argument(self, *args: str, **kwargs: Any) -> Action:
            # JCVI 1.6.6 normalizes help at parse time, after Python 3.14 checks it.
            if isinstance(kwargs.get("help"), str):
                kwargs["help"] = kwargs["help"].replace("%default", "%(default)s")
            return super().add_argument(*args, **kwargs)

    original = quota.OptionParser
    quota.OptionParser = QuotaOptionParser
    try:
        catalog.main()
    finally:
        quota.OptionParser = original


if __name__ == "__main__":
    main()
