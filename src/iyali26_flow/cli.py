"""Command-line interface for the gated two-stage experiment."""

from __future__ import annotations

import argparse
import json
import signal
import sys
from pathlib import Path
from types import FrameType
from typing import Sequence

from .config import load_experiment_config
from .experiment import (
    ExperimentRunner,
    HardTimeoutExceeded,
)
from .fmpe import train_fmpe


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="iyali26-flow",
        description="Gated R4/R1846 parameter-inference experiments for iYali26",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("phase1", "analyze", "train-fmpe"):
        command = subparsers.add_parser(name)
        command.add_argument("--config", type=Path, required=True)
        command.add_argument("--output", type=Path, required=True)
        command.add_argument(
            "--research-root",
            type=Path,
            help="External research workspace; overrides IYALI26_RESEARCH_ROOT",
        )
        if name == "phase1":
            command.add_argument(
                "--no-resume",
                action="store_true",
                help="refuse an existing output instead of resuming it",
            )
    return parser


def _install_alarm(seconds: int):
    if not hasattr(signal, "SIGALRM"):
        return None

    def handler(_signum: int, _frame: FrameType | None) -> None:
        raise HardTimeoutExceeded(
            f"Hard timeout reached after {seconds} seconds; checkpoint is resumable"
        )

    previous = signal.signal(signal.SIGALRM, handler)
    signal.alarm(seconds)
    return previous


def _clear_alarm(previous) -> None:
    if previous is not None:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        config = load_experiment_config(
            args.config,
            repo_root=Path.cwd(),
            research_root=args.research_root,
        )
        previous_alarm = _install_alarm(config.hard_timeout_seconds)
        try:
            if args.command == "phase1":
                runner = ExperimentRunner(
                    config,
                    args.output,
                    resume=not args.no_resume,
                )
                outcome = runner.run_phase1()
                payload = outcome.to_dict()
                exit_code = 124 if outcome.status == "partial_timeout" else 0
            elif args.command == "analyze":
                runner = ExperimentRunner(config, args.output, resume=True)
                outcome = runner.analyze_only()
                payload = outcome.to_dict()
                exit_code = 0
            else:
                payload = train_fmpe(config, args.output)
                exit_code = 0 if payload.get("acceptance_pass") else 1
        finally:
            _clear_alarm(previous_alarm)
    except Exception as exc:
        print(
            json.dumps(
                {"status": "error", "error": str(exc)},
                ensure_ascii=False,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 1
    print(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True))
    return exit_code
