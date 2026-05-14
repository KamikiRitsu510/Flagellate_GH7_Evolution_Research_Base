"""Shared utilities: email config, rate limiting, retry logic."""

import time
import functools
import logging
from typing import Callable, TypeVar

from Bio import Entrez

logger = logging.getLogger(__name__)

# ── Email configuration ──────────────────────────────────────────────────────

_email_set = False


def set_email(email: str) -> None:
    """Set the email address sent to NCBI with every Entrez request.

    NCBI requires a valid email so they can contact you if your usage causes
    problems. Call this once at the start of your script.

    Parameters
    ----------
    email : str
        A valid email address (e.g. ``"you@institution.edu"``).

    Examples
    --------
    >>> from phylofetch import set_email
    >>> set_email("you@example.com")
    """
    global _email_set
    Entrez.email = email
    _email_set = True
    logger.debug("Entrez email set to %s", email)


def _ensure_email() -> None:
    if not _email_set and not Entrez.email:
        raise RuntimeError(
            "No NCBI email configured. Call `phylofetch.set_email('you@example.com')` "
            "before making any requests."
        )


# ── Rate limiter ─────────────────────────────────────────────────────────────

class RateLimiter:
    """Callable decorator that enforces a minimum inter-call interval.

    NCBI limits unauthenticated users to ≤3 requests/second. With an API key
    (set ``Entrez.api_key``) the limit rises to 10/second.

    Parameters
    ----------
    min_interval : float
        Minimum seconds between successive calls. Default ``0.4`` (≈ 2.5 req/s).

    Examples
    --------
    >>> limiter = RateLimiter(0.4)
    >>> limiter.wait()   # blocks until the interval has elapsed
    """

    def __init__(self, min_interval: float = 0.4):
        self.min_interval = min_interval
        self._last_call: float = 0.0

    def wait(self) -> None:
        """Block until at least ``min_interval`` seconds since the last call."""
        elapsed = time.monotonic() - self._last_call
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self._last_call = time.monotonic()

    def __call__(self, fn: Callable) -> Callable:
        """Use as a decorator: ``@limiter``."""
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            self.wait()
            return fn(*args, **kwargs)
        return wrapper


# Shared default limiter (safe for NCBI without API key)
_default_limiter = RateLimiter(min_interval=0.4)


# ── Retry logic ───────────────────────────────────────────────────────────────

F = TypeVar("F", bound=Callable)


def with_retry(fn: F, retries: int = 3, backoff: float = 2.0) -> F:
    """Wrap *fn* so that transient network errors are retried automatically.

    Parameters
    ----------
    fn : callable
        The function to wrap.
    retries : int
        Maximum number of retry attempts after the first failure. Default 3.
    backoff : float
        Multiplicative back-off factor between attempts. Default 2.0.

    Returns
    -------
    callable
        Wrapped function with the same signature as *fn*.
    """
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        delay = 1.0
        last_exc = None
        for attempt in range(retries + 1):
            try:
                return fn(*args, **kwargs)
            except Exception as exc:
                last_exc = exc
                if attempt < retries:
                    logger.warning(
                        "%s failed (attempt %d/%d): %s — retrying in %.1fs",
                        fn.__name__, attempt + 1, retries + 1, exc, delay
                    )
                    time.sleep(delay)
                    delay *= backoff
        raise last_exc  # type: ignore[misc]
    return wrapper  # type: ignore[return-value]
