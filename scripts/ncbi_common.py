"""Minimal helper for calling NCBI E-utilities.

NCBI allows 3 requests per second without an API key and 10 with one.
Set the environment variable to use a key:
    export NCBI_API_KEY=xxxxxxxx
"""

import os
import sys
import threading
import time

import requests
from requests.adapters import HTTPAdapter

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "16s_rRNA_gene_conservation"

API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
MIN_INTERVAL = 0.11 if API_KEY else 0.35  # seconds between requests


class RateLimiter:
    def __init__(self, interval=MIN_INTERVAL):
        self.interval = interval
        self.lock = threading.Lock()
        self.last = 0.0

    def wait(self):
        with self.lock:
            gap = time.time() - self.last
            if gap < self.interval:
                time.sleep(self.interval - gap)
            self.last = time.time()


LIMITER = RateLimiter()


def make_session():
    s = requests.Session()
    adapter = HTTPAdapter(pool_connections=2, pool_maxsize=4, max_retries=0)
    s.mount("https://", adapter)
    s.headers.update({"User-Agent": TOOL + " (python-requests)"})
    return s


def eutils(session, endpoint, params, email=None, attempts=5, timeout=180):
    """Call one E-utilities endpoint via POST, retrying on transient failures."""
    payload = dict(params)
    payload["tool"] = TOOL
    if API_KEY:
        payload["api_key"] = API_KEY
    if email:
        payload["email"] = email

    last = ""
    for i in range(attempts):
        LIMITER.wait()
        try:
            r = session.post("%s/%s.fcgi" % (EUTILS, endpoint), data=payload, timeout=timeout)
        except requests.RequestException as exc:
            last = "%s: %s" % (type(exc).__name__, exc)
        else:
            if r.status_code == 200:
                return r.text
            last = "HTTP %d" % r.status_code
            if r.status_code not in (429, 500, 502, 503, 504):
                break
        wait = min(60, 3 * 2 ** i)
        print("    retrying in %ds (%s)" % (wait, last), file=sys.stderr, flush=True)
        time.sleep(wait)
    raise RuntimeError("E-utilities call failed: %s (%s)" % (endpoint, last))
