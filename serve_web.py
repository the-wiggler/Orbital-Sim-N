#!/usr/bin/env python3
"""
Simple local server for the Emscripten build that:
- Serves the build output directory containing OrbitSimulation.html
- Adds the required cross-origin isolation headers so pthreads work:
  - Cross-Origin-Opener-Policy: same-origin
  - Cross-Origin-Embedder-Policy: require-corp
- Opens the default browser to OrbitSimulation.html on startup

Usage examples:
  python serve_web.py                    # auto-detects build dir and serves on 127.0.0.1:8000
  python serve_web.py --port 8080        # choose a port
  python serve_web.py --root build/Release  # specify directory that contains OrbitSimulation.html

Note: This is intended for local development/testing.
"""

import argparse
import contextlib
import os
import socket
import sys
import threading
import time
import webbrowser
from functools import partial
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
import mimetypes


DEFAULT_HTML = "OrbitSimulation.html"


def guess_root(candidates=None):
    """Try to find a directory that contains OrbitSimulation.html.
    Returns absolute path or None.
    """
    if candidates is None:
        candidates = [
            ".",
            os.path.join("build"),
            os.path.join("build", "Release"),
            os.path.join("build", "Debug"),
            "cmake-build-debug",
            "cmake-build-release",
            os.path.join("cmake-build-debug", "Release"),
            os.path.join("cmake-build-release", "Release"),
        ]

    for rel in candidates:
        candidate = os.path.abspath(rel)
        html_path = os.path.join(candidate, DEFAULT_HTML)
        if os.path.isfile(html_path):
            return candidate
    return None


class COOPCOEPRequestHandler(SimpleHTTPRequestHandler):
    """HTTP handler that injects COOP/COEP headers and sets correct MIME types."""

    # Add headers required by Emscripten pthreads
    def end_headers(self):
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        # Helpful for static assets loaded from same-origin
        self.send_header("Cross-Origin-Resource-Policy", "same-origin")
        # Disable caches for iterative development
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()

    # Optional nicety: redirect "/" to the app HTML
    def do_GET(self):
        if self.path in ("/", ""):
            self.send_response(302)
            self.send_header("Location", f"/{DEFAULT_HTML}")
            self.end_headers()
            return
        return super().do_GET()


def ensure_mime_types():
    # Ensure correct .wasm MIME type
    mimetypes.add_type("application/wasm", ".wasm", strict=True)
    # Common source maps
    mimetypes.add_type("application/json", ".map", strict=False)


def wait_for_port(host: str, port: int, timeout: float = 3.0) -> bool:
    """Wait briefly until the server socket is ready."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
            sock.settimeout(0.2)
            try:
                if sock.connect_ex((host, port)) == 0:
                    return True
            except OSError:
                pass
        time.sleep(0.05)
    return False


def main(argv=None):
    parser = argparse.ArgumentParser(description="Serve OrbitSimulation.html with COOP/COEP headers")
    parser.add_argument("--host", default="127.0.0.1", help="Interface to bind (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind (default: 8000)")
    parser.add_argument("--root", default=None, help="Directory containing OrbitSimulation.html")
    parser.add_argument("--no-open", action="store_true", help="Do not open a browser automatically")
    args = parser.parse_args(argv)

    # Resolve root directory
    root = args.root or guess_root()
    if not root:
        print(f"[serve_web] Could not find {DEFAULT_HTML}. Specify --root pointing to your build output.", file=sys.stderr)
        return 2

    root = os.path.abspath(root)
    html_path = os.path.join(root, DEFAULT_HTML)
    if not os.path.isfile(html_path):
        print(f"[serve_web] {DEFAULT_HTML} not found in {root}", file=sys.stderr)
        return 2

    ensure_mime_types()

    handler_cls = partial(COOPCOEPRequestHandler, directory=root)
    server = ThreadingHTTPServer((args.host, args.port), handler_cls)

    url = f"http://{args.host}:{args.port}/{DEFAULT_HTML}"
    print(f"[serve_web] Serving {root}")
    print("[serve_web] Headers: COOP=same-origin, COEP=require-corp (threads enabled)")
    print(f"[serve_web] Open: {url}")

    # Start server in background thread so we can open browser after bind
    t = threading.Thread(target=server.serve_forever, name="serve_web", daemon=True)
    t.start()

    opened = False
    if not args.no_open:
        if wait_for_port(args.host, args.port, timeout=3.0):
            try:
                webbrowser.open(url)
                opened = True
            except Exception:
                opened = False

    if not opened and not args.no_open:
        print("[serve_web] Could not auto-open browser. Open the URL above manually.")

    try:
        # Block the main thread while server runs
        while t.is_alive():
            time.sleep(0.2)
    except KeyboardInterrupt:
        print("\n[serve_web] Shutting down...")
    finally:
        with contextlib.suppress(Exception):
            server.shutdown()
            server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
