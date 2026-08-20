"""Tests for pipeline output capture."""

from __future__ import annotations

import time

from rich.console import Console

from lightsuite.reporter import CallbackReporter, capture_pipeline_output


def test_capture_pipeline_output_forwards_rich_console() -> None:
    messages: list[str] = []
    reporter = CallbackReporter(on_message=messages.append)
    with capture_pipeline_output(reporter):
        Console().print("slice progress")
    assert messages.count("slice progress") == 1


def test_capture_pipeline_output_forwards_stdout() -> None:
    messages: list[str] = []
    reporter = CallbackReporter(on_message=messages.append)
    with capture_pipeline_output(reporter):
        print("plain stdout line")
    assert any("plain stdout line" in message for message in messages)


def test_report_step_progress_forwards_to_reporter() -> None:
    messages: list[str] = []
    reporter = CallbackReporter(on_message=messages.append)
    t0 = time.perf_counter()
    with capture_pipeline_output(reporter):
        from lightsuite.reporter import report_step_progress

        report_step_progress(1, 10, label="channel 1", t0=t0, unit="slice")
    assert len(messages) == 1
    assert "channel 1" in messages[0]
    assert "slice 1/10" in messages[0]
    assert "elapsed" in messages[0]
    assert "left" in messages[0]


def test_format_duration() -> None:
    from lightsuite.reporter import format_duration

    assert format_duration(12) == "12s"
    assert format_duration(65) == "1m 05s"
    assert format_duration(3725) == "1h 02m"
