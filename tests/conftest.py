"""Shared fixtures and import-time mocks for MOEA-FIND tests.

objectives.py imports synhydro at module level, so we must stub it out
before any test module causes that import to run.
"""

import sys
from unittest.mock import MagicMock

# Stub synhydro so objectives.py (and any other src module) can be imported
# without SynHydro installed.
_synhydro = MagicMock()
sys.modules.setdefault("synhydro", _synhydro)
sys.modules.setdefault("synhydro.droughts", _synhydro.droughts)
sys.modules.setdefault("synhydro.droughts.ssi", _synhydro.droughts.ssi)
