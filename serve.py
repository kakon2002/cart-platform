"""Start the platform's HTTP server.

A single entry point so the server module is imported once under its real
name. Running `python -m car_pipeline.api.server` executes that file as
__main__ and imports it again when the adapter reaches it, leaving two copies
whose exception classes are different objects -- and a refusal raised by one
is then not caught by the other.
"""

from __future__ import annotations

import argparse
import os

from car_pipeline.api import server


def main() -> int:
    """Parse the binding and serve."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--host", default=os.environ.get("HOST") or "127.0.0.1")
    parser.add_argument("--port", type=int,
                        default=int(os.environ.get("PORT") or 8000))
    options = parser.parse_args()
    server.serve(options.host, options.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
