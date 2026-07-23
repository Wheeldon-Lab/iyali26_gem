"""
http_utils.py — shared HTTP helper with retry and exponential backoff.
"""

import logging
import time

import requests

logger = logging.getLogger(__name__)

_RETRYABLE_STATUS = {429, 500, 502, 503, 504}
_RETRYABLE_EXC   = (requests.exceptions.ConnectionError, requests.exceptions.Timeout)


def _request_with_retry(method: str, url: str, **kwargs) -> requests.Response:
    """
    Execute an HTTP request with up to 3 retries on transient failures.

    Retries on: ConnectionError, Timeout, and HTTP 429/500/502/503/504.
    Backoff: 1 s, 2 s, 4 s between attempts.
    On HTTP 429, the Retry-After header is respected if present.

    Non-retryable failures (e.g. 400, 404) raise immediately after the first attempt.
    """
    delays = [1, 2, 4]
    last_exc: Exception | None = None

    for attempt, delay in enumerate(delays + [None], start=1):
        try:
            resp = requests.request(method, url, **kwargs)
        except _RETRYABLE_EXC as exc:
            last_exc = exc
            if delay is None:
                break
            logger.debug(f"  HTTP {method} {url}: {exc} — retry {attempt}/3 in {delay}s")
            time.sleep(delay)
            continue

        if resp.status_code not in _RETRYABLE_STATUS:
            return resp   # success or non-retryable error — let caller handle

        if delay is None:
            return resp   # exhausted retries, return the last response

        wait = delay
        if resp.status_code == 429:
            try:
                wait = int(resp.headers.get("Retry-After", delay))
            except (TypeError, ValueError):
                pass
        logger.debug(
            f"  HTTP {method} {url}: status {resp.status_code} — "
            f"retry {attempt}/3 in {wait}s"
        )
        time.sleep(wait)

    if last_exc is not None:
        raise last_exc
    return resp  # unreachable, but satisfies type checker
