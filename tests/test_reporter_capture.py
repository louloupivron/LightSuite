"""Tests for pipeline output capture."""

from __future__ import annotations

from rich.console import Console

from lightsuite.reporter import CallbackReporter, capture_pipeline_output


def test_capture_pipeline_output_forwards_rich_console() -> None:
    messages: list[str] = []
    reporter = CallbackReporter(on_message=messages.append)
    with capture_pipeline_output(reporter):
        Console().print("slice progress")
    assert any("slice progress" in message for message in messages)


def test_capture_pipeline_output_forwards_stdout() -> None:
    messages: list[str] = []
    reporter = CallbackReporter(on_message=messages.append)
    with capture_pipeline_output(reporter):
        print("plain stdout line")
    assert any("plain stdout line" in message for message in messages)
