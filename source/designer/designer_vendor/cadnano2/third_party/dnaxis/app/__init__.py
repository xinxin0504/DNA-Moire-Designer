"""DNAxiS package bootstrap.

cadnano uses the routing engine headlessly and deliberately does not import
the Flask web interface.  The original behaviour remains available when the
environment switch is absent.
"""

import os

if os.environ.get("DNAXIS_HEADLESS") == "1":
    app = None
    from app import config
else:
    from flask import Flask
    app = Flask(__name__)
    from app import config
    from app import views
    from app import forms
